#!/usr/bin/env python3
"""Compare source-basis XML HAR to the matching library path on supplied data.

No calculation data are downloaded, changed in place, or checked into Git.
The user supplies the XML, explicit source basis, CIF, HKL and library name.
The checker runs otherwise identical inputs in isolated fresh directories.
"""
from __future__ import annotations
import argparse
from decimal import Decimal, InvalidOperation
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys


def equivalent_phase_rows(original: bytes, imported: bytes) -> int:
    """Allow exact whole-turn phase representations, never a numeric tolerance.

    All other rows/fields must agree. Retain byte hashes and report the number
    of equivalent phase rows instead of describing such files as identical.
    """
    left_lines, right_lines = original.decode().splitlines(), imported.decode().splitlines()
    if len(left_lines) != len(right_lines):
        raise ValueError("Reflection outputs have different line counts")
    headers, phase_index, changed = [], None, 0
    for left, right in zip(left_lines, right_lines):
        stripped = left.strip()
        control = stripped.lower()
        # A new block/frame or explicit loop terminator ends this column map.
        # Never reuse a previous reflection loop's phase tag on later rows.
        if (control in ("loop_", "stop_")
                or control.startswith(("data_", "save_"))):
            headers, phase_index = [], None
        elif stripped.startswith("_"):
            headers.append(stripped)
            if stripped == "_refln_phase_calc":
                phase_index = len(headers) - 1
        if left == right:
            continue
        a, b = left.split(), right.split()
        if (phase_index is None or len(headers) != len(set(headers))
                or len(a) != len(headers) or len(b) != len(headers)
                or not all(name in headers for name in
                           ("_refln_index_h", "_refln_index_k", "_refln_index_l"))):
            raise ValueError("Difference outside a tagged reflection phase row")
        if any(x != y for i, (x, y) in enumerate(zip(a, b)) if i != phase_index):
            raise ValueError("Non-phase reflection fields differ")
        for name in ("_refln_index_h", "_refln_index_k", "_refln_index_l"):
            if not re.fullmatch(r"[-+]?\d+", a[headers.index(name)]):
                raise ValueError("Invalid reflection index")
        try:
            p, q = Decimal(a[phase_index]), Decimal(b[phase_index])
            if not p.is_finite() or not q.is_finite() or (Fraction(p) - Fraction(q)) % 360 != 0:
                raise ValueError("Phases are not equal modulo 360 degrees")
        except InvalidOperation as error:
            raise ValueError("Invalid reflection phase") from error
        changed += 1
    if not changed:
        raise ValueError("Nonidentical bytes without an equivalent phase difference")
    return changed


def validate_existing(root: Path) -> dict[str, object]:
    """Read-only, repeatable validation of an already completed isolated pair."""
    texts = {}
    for mode in ("library", "companion"):
        output = root / mode / "stdout"
        text = output.read_text(errors="replace")
        if "Structure fit converged." not in text:
            raise RuntimeError(f"{mode} did not converge: {output}")
        if re.search(r"\b(?:NaN|Inf|Infinity)\b", text, re.I):
            raise RuntimeError(f"{mode} reported a nonfinite value: {output}")
        if ("Initial density matrix for refinement" in text
                or re.search(r"(?m)^SCF calculation\s*$", text)):
            raise RuntimeError(f"{mode} started a molecular SCF for imported periodic density")
        texts[mode] = text
    final = {key: value.rsplit("Structure fit converged.", 1)[1].split("Wall-clock time", 1)[0]
             for key, value in texts.items()}
    numbers = {key: [float(x.replace("D", "E")) for x in
                    re.findall(r"(?<![A-Za-z_])[-+]?\d*\.\d+(?:[EeDd][-+]?\d+)?", value)]
               for key, value in final.items()}
    left, right = numbers["library"], numbers["companion"]
    if len(left) != len(right) or not left:
        raise RuntimeError("Final reports have inconsistent numeric fields")
    maximum = max(abs(a-b) for a, b in zip(left, right))
    # The unchanged printed-report tolerance is two final units of 1e-6.
    if maximum > 2e-6:
        raise RuntimeError(f"Source/library report differs by {maximum:.12g}; evidence {root}")
    artifacts, artifact_hashes, phase_rows = {}, {}, {}
    for suffix in ("archive.fcf", "archive.fco", "fcf6"):
        name = "xml_companion_real." + suffix
        left_file, right_file = root / "library" / name, root / "companion" / name
        if not left_file.is_file() or not right_file.is_file():
            raise RuntimeError(f"Missing calculated/observed structure-factor output {name}")
        original, imported = left_file.read_bytes(), right_file.read_bytes()
        artifacts[name] = original == imported
        artifact_hashes[name] = {"library": hashlib.sha256(original).hexdigest(),
                                 "companion": hashlib.sha256(imported).hexdigest()}
        if not artifacts[name]:
            try:
                if suffix == "archive.fco":
                    raise ValueError("FCO outputs must be byte-identical")
                phase_rows[name] = equivalent_phase_rows(original, imported)
            except ValueError as error:
                raise RuntimeError(f"Source/library reflection output differs: {name}; "
                                   f"{error}; evidence {root}") from error
        else:
            phase_rows[name] = 0
    return {"numeric_fields": len(left), "maximum_report_difference": maximum,
            "molecular_scf_detected": False,
            "reflection_outputs_byte_identical": artifacts,
            "reflection_outputs_equivalent": True,
            "equivalent_phase_representation_rows": phase_rows,
            "reflection_output_sha256": artifact_hashes, "passed": True}


def main() -> int:
    if "--check-existing" in sys.argv:
        parser = argparse.ArgumentParser(description="Read-only recheck of a completed source/library pair")
        parser.add_argument("--check-existing", required=True, type=Path)
        args = parser.parse_args()
        print(json.dumps(validate_existing(args.check_existing.resolve()), indent=2))
        return 0
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--program", required=True, type=Path)
    parser.add_argument("--xml", required=True, type=Path)
    parser.add_argument("--source-basis", required=True, type=Path)
    parser.add_argument("--cif", required=True, type=Path)
    parser.add_argument("--hkl", required=True, type=Path)
    parser.add_argument("--basis-directory", required=True, type=Path)
    parser.add_argument("--basis-name", required=True)
    parser.add_argument("--space-group", required=True)
    parser.add_argument("--stockholder", choices=("cluster", "periodic"), default="cluster")
    parser.add_argument("--work-directory", required=True, type=Path)
    parser.add_argument("--timeout", type=int, default=900,
                        help="Seconds per calculation (default: 900)")
    args = parser.parse_args()
    if args.timeout <= 0:
        parser.error("--timeout must be positive")
    root = args.work_directory.resolve()
    root.mkdir(parents=True, exist_ok=True)
    if any(root.iterdir()):
        raise RuntimeError(f"Refusing to overwrite nonempty directory {root}")
    paths = {name: getattr(args, name).resolve() for name in ("program", "xml", "source_basis", "cif", "hkl", "basis_directory")}
    manifest = {}
    for name in ("xml", "source_basis", "cif", "hkl"):
        path = paths[name]
        with path.open("rb") as stream:
            checksum = hashlib.sha256()
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                checksum.update(chunk)
            digest = checksum.hexdigest()
        manifest[name] = {"path": str(path), "sha256": digest}
    environment = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    for mode in ("library", "companion"):
        work = root / mode
        work.mkdir()
        os.symlink(paths["xml"], work / "input.xml")
        shutil.copy2(paths["cif"], work / "input.cif")
        shutil.copy2(paths["hkl"], work / "input.hkl")
        shutil.copy2(paths["source_basis"], work / "source.basis")
        basis = (f"basis_name= {args.basis_name}"
                 if mode == "library" else "c23_xml_basis_file= source.basis")
        entry = f'''{{
   ! Independent old-library/source-companion paths, otherwise identical.
   name= xml_companion_real
   ! Expose an accidental molecular-SCF fallback; do not execute an SCF.
   scf_data= {{ kind= rhf output= true }}
   CIF= {{ file_name= input.cif }}
   ! Thakkar is the partition proatom, not the imported orbital basis.
   basis_directory= {paths['basis_directory']}
   slaterbasis_name= Thakkar
   {basis}
   charge= 0
   multiplicity= 1
   c23_XML_file_name= input.xml
   process_cif_and_c23_xml
   crystal= {{
      spacegroup= {{ hermann_mauguin_symbol= "{args.space_group}" }}
      xray_data= {{
         partition_model= oc-crystal23
         stockholder_model= {args.stockholder}
         optimise_extinction= false
         correct_dispersion= no
         optimise_scale_factor= true
         wavelength= 0.71073 Angstrom
         REDIRECT input.hkl
         f_sigma_cutoff= 4
         tol_for_shift_on_esd= 0.01
         refine_H_U_iso= no
         max_iterations= 30
         show_fit_output= false
         show_fit_results= true
      }}
   }}
   phar_defragment
   ha_fit
}}
'''
        (work / "stdin").write_text(entry)
        with (work / "driver.log").open("w") as log:
            run = subprocess.run([str(paths["program"])], cwd=work, env=environment,
                                 stdout=log, stderr=subprocess.STDOUT, timeout=args.timeout)
        text = (work / "stdout").read_text(errors="replace")
        if run.returncode or "Structure fit converged." not in text:
            raise RuntimeError(f"{mode} did not converge: {work / 'stdout'}")
        if re.search(r"\b(?:NaN|Infinity)\b", text, re.I):
            raise RuntimeError(f"{mode} reported a nonfinite value: {work / 'stdout'}")
        print(f"PASS {mode} converged", flush=True)
    manifest.update(validate_existing(root))
    manifest["stockholder"] = args.stockholder
    (root / "comparison.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"PASS real XML source/library comparison: {manifest['numeric_fields']} numbers, "
          f"maxdelta={manifest['maximum_report_difference']:g}; {root}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (RuntimeError, subprocess.TimeoutExpired) as error:
        print(f"FAIL: {error}", file=sys.stderr)
        raise SystemExit(1)
