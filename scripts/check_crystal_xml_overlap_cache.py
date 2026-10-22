#!/usr/bin/env python3
"""Exercise the real XML-import validator after a cache has been populated.

There are no refreshed numerical references: unchanged inputs must avoid the
integral calculation, changed valid inputs must repeat it, and invalid inputs
must fail before a cache hit can conceal them. The analytic source fixture is
copied, never edited in place.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import time


VALID = (
    "unchanged", "reread", "geometry", "element", "source_overlap",
    "order_flag", "global_order_flag", "spherical_flag", "shell_update", "exponent", "contraction",
)
INVALID = {
    "nan_overlap": r"CRYSTAL.*nonfinite.*overlap",
    "inf_overlap": r"CRYSTAL.*nonfinite.*overlap",
    "bad_dimension": r"CRYSTAL.*overlap dimension",
    "nan_position": r"CRYSTAL XML.*nonfinite.*(?:basis|geometry)",
    "nan_exponent": r"CRYSTAL XML.*nonfinite.*(?:basis|geometry)",
    "inf_contraction": r"CRYSTAL XML.*nonfinite.*(?:basis|geometry)",
    "bad_atom_map": r"CRYSTAL XML.*(?:atom/shell|atom.shell).*mapping",
    "bad_ao_map": r"CRYSTAL XML.*AO shell ranges",
    "bad_basis_map": r"CRYSTAL XML.*basis-shell mapping",
    "bad_overlap": r"CRYSTAL.*(?:basis overlap validation failed|basis differs)",
}
FIXTURE_FILES = ("stdin", "input.cif", "input.xml", "source.basis")


def exercise(program: Path, fixture: Path, work: Path) -> list[dict]:
    results = []
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    for scenario in (*VALID, *INVALID):
        folder = work / scenario
        folder.mkdir(parents=True, exist_ok=False)
        for name in FIXTURE_FILES:
            shutil.copy2(fixture / name, folder / name)
        started = time.perf_counter()
        run = subprocess.run(
            [str(program), scenario], cwd=folder, env=env,
            capture_output=True, text=True, timeout=60,
        )
        output_file = folder / "stdout"
        output = (output_file.read_text(errors="replace") if output_file.exists() else "")
        output += "\n" + run.stdout + "\n" + run.stderr
        if re.search(r"SIGSEGV|Segmentation fault|invalid memory reference", output, re.I):
            raise AssertionError(f"{scenario}: crashed instead of handling the cached state\n{output}")
        marker = f"CACHE MUTATION: {scenario}"
        if marker not in output:
            raise AssertionError(f"{scenario}: cache was not populated before mutation\n{output}")
        after_mutation = output.split(marker, 1)[1]
        accepted = f"CACHE SCENARIO ACCEPTED: {scenario}" in after_mutation
        if scenario in INVALID:
            if accepted or re.search(INVALID[scenario], after_mutation, re.IGNORECASE) is None:
                raise AssertionError(
                    f"{scenario}: expected production rejection, not a hit or unrelated crash\n{output}"
                )
        elif run.returncode != 0 or not accepted:
            raise AssertionError(f"{scenario}: cache/invalidation assertion failed\n{output}")
        results.append({
            "scenario": scenario,
            "outcome": "rejected" if scenario in INVALID else "accepted",
            "returncode": run.returncode,
            "elapsed_seconds": time.perf_counter() - started,
        })
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--program", type=Path, required=True)
    parser.add_argument("--fixture", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, help="Optional new directory retaining every stdout")
    args = parser.parse_args()
    program, fixture = args.program.resolve(), args.fixture.resolve()
    if args.work_dir:
        work = args.work_dir.resolve()
        work.mkdir(parents=True, exist_ok=False)
        results = exercise(program, fixture, work)
        (work / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    else:
        with tempfile.TemporaryDirectory(prefix="tonto-xml-overlap-cache-") as temporary:
            results = exercise(program, fixture, Path(temporary))
    print(f"XML overlap cache: {len(results)} production integration cases passed")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
