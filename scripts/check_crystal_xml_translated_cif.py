#!/usr/bin/env python3
"""Check periodic CIF-to-XML matching without changing the source wavefunction.

The analytic H/He/H fixture has distinct bases on its two H sites. Only CIF
representatives are translated. The source XML and companion basis stay fixed.
"""
import argparse
import json
from pathlib import Path
import tempfile

from check_crystal_xml_companion_basis import demand, run_case


def shifted_cif(text, shifts, displaced=False):
    lines, count = [], 0
    labels = ("H1", "He1", "H2")
    for line in text.splitlines():
        fields = line.split()
        if len(fields) == 8 and fields[0] in labels:
            index = labels.index(fields[0])
            position = [float(x)+t for x,t in zip(fields[2:5],shifts[index])]
            if displaced and index == 0:
                position[0] += .02  # 0.2 bohr, outside the unchanged 0.1-bohr tolerance
            fields[2:5] = [f"{x:.12f}" for x in position]
            line = " ".join(fields)
            count += 1
        lines.append(line)
    demand(count == 3, "Expected exactly the three analytic CIF sites")
    return "\n".join(lines)+"\n"


def check(program, fixture, root):
    original = {name:(fixture/name).read_text(encoding="ascii")
                for name in ("stdin","input.cif","input.xml","source.basis")}
    cases = {
        "untranslated": ((0,0,0),(0,0,0),(0,0,0)),
        "negative_and_positive_cells": ((-4,3,-2),(5,-6,4),(-7,-8,9)),
        "opposite_cells": ((4,-3,2),(-5,6,-4),(7,8,-9)),
    }
    results = []
    for label, shifts in cases.items():
        files = dict(original)
        files["input.cif"] = shifted_cif(original["input.cif"],shifts)
        results.append(run_case(program,root,label,files))
        print("PASS",label,flush=True)
    files = dict(original)
    files["input.cif"] = shifted_cif(original["input.cif"],
                                   cases["negative_and_positive_cells"],displaced=True)
    results.append(run_case(program,root,"translated_but_wrong_site",files,
                            r"A cif atom is not found among.*XML atoms"))
    print("PASS translated_but_wrong_site rejected",flush=True)
    for result in results:
        folder = Path(result["path"])
        for name in ("input.xml","source.basis"):
            demand((folder/name).read_text() == original[name],f"Changed {name}")
        demand(not list(folder.glob("*.GRED")), "XML test must not use a GRED file")
    (root/"results.json").write_text(json.dumps(results,indent=2)+"\n")
    print("All 4 CIF/XML translation checks passed")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--program",type=Path,required=True)
    parser.add_argument("--fixture",type=Path,required=True)
    parser.add_argument("--work-directory",type=Path)
    args = parser.parse_args()
    if args.work_directory:
        root = args.work_directory.resolve()
        root.mkdir(parents=True,exist_ok=False)
        check(args.program.resolve(),args.fixture.resolve(),root)
    else:
        with tempfile.TemporaryDirectory(prefix="tonto-xml-translation-") as folder:
            check(args.program.resolve(),args.fixture.resolve(),Path(folder))


if __name__ == "__main__":
    main()
