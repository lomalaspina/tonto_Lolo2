#!/usr/bin/env python3
"""Exercise the XML companion-basis import without an external calculation.

The committed fixture is analytic, not a CRYSTAL-produced wavefunction.  This
checker tests parsing, atom-specific basis assignment and fail-closed overlap
validation.  It never rewrites the committed fixture or a golden stdout.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile


def demand(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def run_case(program: Path, root: Path, label: str, files: dict[str, str],
             failure: str | None = None) -> dict[str, object]:
    work = root / label
    work.mkdir()
    for name, content in files.items():
        (work / name).write_text(content, encoding="ascii")
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    with (work / "driver.log").open("w") as log:
        process = subprocess.run([str(program)], cwd=work, env=env,
                                 stdout=log, stderr=subprocess.STDOUT, timeout=120)
    output = (work / "stdout").read_text(errors="replace") if (work / "stdout").exists() else ""
    driver = (work / "driver.log").read_text(errors="replace")
    combined = output + "\n" + driver
    demand(not re.search(r"SIGSEGV|Segmentation fault|invalid memory reference", combined),
           f"{label}: crashed rather than handling the input")
    if failure:
        demand(process.returncode != 0, f"{label}: invalid input was accepted")
        demand(re.search(failure, combined, re.I) is not None,
               f"{label}: expected {failure!r}; inspect {work / 'stdout'}")
    else:
        demand(process.returncode == 0, f"{label}: failed; inspect {work / 'stdout'}")
        match = re.search(r"Maximum CRYSTAL23/Tonto basis overlap difference[ .=]*([-+0-9.EeDd]+)", output)
        demand(match is not None, f"{label}: central-overlap validation did not run")
        error = float(match.group(1).replace("D", "E").replace("d", "e"))
        demand(math.isfinite(error) and error < 1e-10,
               f"{label}: analytic central overlap differs by {error}")
    return {"case": label, "returncode": process.returncode,
            "expected": "rejection" if failure else "success", "path": str(work)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--program", required=True, type=Path)
    parser.add_argument("--fixture", required=True, type=Path)
    parser.add_argument("--work-directory", type=Path)
    args = parser.parse_args()
    program, fixture = args.program.resolve(), args.fixture.resolve()
    if args.work_directory:
        root = args.work_directory.resolve()
        root.mkdir(parents=True, exist_ok=True)
        demand(not any(root.iterdir()), f"Refusing to overwrite nonempty work directory: {root}")
    else:
        root = Path(tempfile.mkdtemp(prefix="xml-companion-basis-"))
    files = {name: (fixture / name).read_text(encoding="ascii")
             for name in ("stdin", "input.cif", "input.xml", "source.basis")}
    results = []

    def case(label: str, changes: dict[str, str] | None = None,
             failure: str | None = None) -> None:
        data = dict(files)
        data.update(changes or {})
        results.append(run_case(program, root, label, data, failure))
        print(f"PASS {label}", flush=True)

    case("same_element_distinct_bases")
    # Only the named XML and explicit companion supply source data. Unrelated
    # native exports in a calculation directory must never become fallbacks.
    case("unrelated_native_exports_ignored", {
        name: "INVALID DECOY: not input to the XML interface\n"
        for name in ("GRED.DAT", "KRED.DAT", "GenerateXML_dat.GRED", "GenerateXML_dat.KRED")
    })
    no_fock = re.sub(r"<DIRECT_FOCK_KOHN-SHAM_MATRIX>.*?</DIRECT_FOCK_KOHN-SHAM_MATRIX>", "",
                     files["input.xml"], flags=re.S)
    demand(no_fock != files["input.xml"], "Fixture has no optional Fock block")
    reread = files["stdin"].rstrip()[:-1] + """
   ! A second import must discard optional matrices from the first file.
   c23_XML_file_name= second.xml
   process_cif_and_c23_xml
}
"""
    case("reread_without_optional_Fock", {"stdin": reread, "second.xml": no_fock})
    output = (root / "reread_without_optional_Fock" / "stdout").read_text()
    availability = re.findall(r"CRYSTAL XML companion Fock matrix available[ .=]*([TF])", output)
    demand(availability == ["T", "F"],
           f"Second import reused stale Fock information: {availability}")
    demand(output.count("Maximum CRYSTAL23/Tonto basis overlap difference") == 2,
           "Both imports must validate their independently imported overlap")
    # Reorder companion blocks without changing ATOM_MAP: source atom order
    # cannot be inferred from element or from the order of definition blocks.
    block1 = "1 1\n0 0 1 1.0 1.0\n0.7 1.0\n"
    block2 = "101 1\n0 0 1 1.0 1.0\n1.3 1.0\n"
    swapped = files["source.basis"].replace(block1 + block2, block2 + block1)
    demand(swapped != files["source.basis"], "Fixture blocks changed; update explicit reorder control")
    case("reordered_basis_definitions", {"source.basis": swapped})
    wrong_map = files["source.basis"].replace("1 1\n2 2\n3 101\n", "1 101\n2 2\n3 1\n")
    demand(wrong_map != files["source.basis"], "Source-atom map changed; update swap control")
    case("wrong_same_element_atom_map", {"source.basis": wrong_map},
         r"basis differs.*overlap|overlap.*basis")
    case("same_AO_wrong_exponent", {"source.basis": files["source.basis"].replace("0.7 1.0", "0.9 1.0")},
         r"basis differs.*overlap|overlap.*basis")
    # Exact primitive-Gaussian control: make He-H2 overlap 2.5e-4, then
    # compare its accurate source value with a deliberately screened zero.
    # Do not accept a source zero merely because the missing overlap is small.
    alpha, beta = 1.1, 1.3
    prefactor = (4 * alpha * beta / (alpha + beta)**2)**0.75
    distance2 = -math.log(2.5e-4 / prefactor) * (alpha + beta) / (alpha * beta)
    h2_x = 2 + math.sqrt(distance2 - 4)
    s_h1_h2 = (4 * 0.7 * beta / (0.7 + beta)**2)**0.75 * math.exp(
        -0.7 * beta / (0.7 + beta) * (h2_x - 2)**2)
    accurate = files["input.xml"].replace("4.0 2.0 2.0", f"{h2_x:.17g} 2.0 2.0")
    accurate = accurate.replace("0.15096111232544537 0.0084643405576893919",
                                f"{s_h1_h2:.17g} 0.00025")
    small_pair_cif = files["input.cif"].replace("H2 H 0.4", f"H2 H {h2_x / 10:.17g}")
    case("small_overlap_accurate_source", {"input.xml": accurate, "input.cif": small_pair_cif})
    screened = accurate.replace("0.00025", "0.0")
    case("small_overlap_screened_source", {"input.xml": screened, "input.cif": small_pair_cif},
         r"source.zero.*overlap|source screening")
    diagnostic = (root / "small_overlap_screened_source" / "stdout").read_text()
    demand("confined to source-zero entries" in diagnostic and "tighter TOLINTEG" in diagnostic,
           "Screened control did not explain the source-accuracy problem")
    demand(re.search(r"Source-zero entries above tolerance[ .]*2", diagnostic) is not None,
           "Screened control did not count the symmetric pair")
    case("small_overlap_wrong_basis", {"input.xml": accurate, "input.cif": small_pair_cif,
         "source.basis": files["source.basis"].replace("0.7 1.0", "0.9 1.0")},
         r"basis differs.*overlap|overlap.*basis")
    case("wrong_AO_count", {"source.basis": files["source.basis"].replace(block1,
         "1 2\n0 0 1 1.0 1.0\n0.7 1.0\n0 0 1 0.0 1.0\n0.3 1.0\n")},
         r"AO|basis count|basis.*dimension")
    no_overlap = re.sub(r"<DIRECT_OVERLAP_MATRIX>.*?</DIRECT_OVERLAP_MATRIX>", "",
                        files["input.xml"], flags=re.S)
    demand(no_overlap != files["input.xml"], "Fixture has no central-overlap block")
    case("missing_overlap", {"input.xml": no_overlap}, r"overlap")
    case("reread_missing_overlap", {"stdin": reread, "second.xml": no_overlap}, r"overlap")
    case("missing_atom_map", {"source.basis": files["source.basis"].split("ATOM_MAP")[0]},
         r"ATOM_MAP|ambiguous|alias")
    case("missing_atom_map_row", {"source.basis": files["source.basis"].replace("3 101\n", "")},
         r"ATOM_MAP|atom.map|end|EOF")
    case("invalid_exponent", {"source.basis": files["source.basis"].replace("0.7 1.0", "-0.7 1.0")},
         r"exponent|positive")
    case("normalization_overflow", {"source.basis": files["source.basis"].replace("0.7 1.0", "1e300 1.0")},
         r"contraction|normaliz|finite|overlap")
    case("unsupported_scaled_shell", {"source.basis": files["source.basis"].replace("0 0 1 1.0 1.0", "0 0 1 1.0 1.2", 1)},
         r"scal|1\.0|one")
    case("unknown_AO_label", {"input.xml": files["input.xml"].replace("\n S\n", "\n PX\n", 1)},
         r"AO|label|shell|ordering")
    for matrix, tag in (("density", "DENSITY"), ("overlap", "OVERLAP"),
                        ("Fock", "FOCK_KOHN-SHAM")):
        for bad_number in ("NaN", "Infinity"):
            pattern = rf"(<DIRECT_{tag}_MATRIX__IVDL\.1[^>]*>\s*)[-+0-9.EeDd]+"
            bad_xml, count = re.subn(pattern, lambda found: found.group(1) + bad_number,
                                    files["input.xml"], count=1)
            demand(count == 1, f"Missing XML matrix {matrix}")
            case(f"nonfinite_{matrix}_{bad_number}", {"input.xml": bad_xml},
                 r"finite|number|invalid|NaN|Infinity")
    (root / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    print(f"All {len(results)} XML companion-basis checks passed; evidence: {root}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (AssertionError, subprocess.TimeoutExpired) as error:
        print(f"FAIL: {error}", file=sys.stderr)
        raise SystemExit(1)
