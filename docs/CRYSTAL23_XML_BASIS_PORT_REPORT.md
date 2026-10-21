# CRYSTAL23 XML companion-basis work record

## Scope and current state

21 October 2026. Companion import and source-screening controls pass locally.
The guarded full-XML L-alanine HAR pairs have completed, as has a new Urea POB
source-screening control and fixed-density HAR pair using the exact received
experiment. Historical incomplete attempts are retained below; the final
sections record their completed follow-up and GRED-free publication checks.
Repository history records publication; unresolved scientific cases remain
explicitly separate from the tested implementation.

The new reader is intended to take an explicit basis definition from the same
CRYSTAL calculation as the XML matrices. It must not infer missing exponents or
contractions from the overlap matrix. The overlap comparison remains a check.

This checkout is the corrected public XML code on branch
`xml-periodic-stockholder-master`, initially at `c81eab01`, with `origin/master`
as its upstream. Its preceding XML-correction commit is `d4779aba`. Only XML
companion-basis changes and their tests and documentation belong in this port.

## Baseline

The checkout had no tracked or untracked changes before this work. The existing
`build-focused` tree uses gfortran-14 and `RELEASE`, with no dead-code purge.
Here the executable target is `tonto` and the library target is `tonto_lib`.

Command:

```sh
ctest --test-dir build-focused -L short --output-on-failure -j2
```

Result: **70/70 short tests passed**, in 19.04 seconds. This is the unchanged
baseline, not evidence that the new reader has been implemented or tested.

## Changes

- Added `native/crystal_xml_basis_io.f90` and its independent parser tests.
- Added the optional `c23_xml_basis_file` keyword and an XML-only basis import.
- Kept one basis per source atom and retained the existing matrix permutation.
- Added AO-label checks and required an early central-overlap check.
- Kept the old library path and its output unchanged when the option is absent.
- Added guards against invalid, overflowing or underflowing contractions.
- Added finite-number checks for imported density, overlap and Fock matrices
  and the overlap constructed from the imported basis.
- Checked that a second XML import cannot reuse an absent overlap or Fock
  matrix from the first import.
- Added the analytic accepted/rejected-input fixture and semantic checker.
- Added a real-data XML library/companion comparison script. No experimental
  input files or large XML files were added to this repository.
- Added user documentation linked from `docs/DOCUMENTATION.md`.

## Final short-test validation

The finalized `tonto` build succeeded. The executable SHA-256 is
`21fb0d71009f80c516011ee312ab650509e855d323f66a01c68dfbd34c2bfb7e`.

Commands used:

```sh
cmake -S . -B build-focused
cmake --build build-focused --target tonto run_crystal_xml_basis_io -j2
python3 scripts/check_crystal_xml_companion_basis.py \
  --program build-focused/tonto \
  --fixture tests/short/crystal_xml_companion_basis \
  --work-directory build-focused/xml-companion-contract-final-20261018
build-focused/run_crystal_xml_basis_io
ctest --test-dir build-focused -L short --output-on-failure -j2
```

- **75/75** native parser checks pass.
- **20/20** semantic import checks pass. These cover two different bases on
  the same element, atom/AO permutation, reordered definitions, wrong source
  atom maps, same-size wrong bases, unsupported inputs, numeric overflow,
  non-finite matrix entries and repeated XML imports.
- The independent analytic three-atom overlap differs by at most
  **2.2204460492503131e-16**. The diagnostic prints enough digits to verify
  the test's `1e-10` tolerance.
- The new normal-design short-test `stdout` was installed only after these
  checks passed. It is a successful import, not a failed output.
- **73/73 short tests pass**, in 19.25 seconds. No existing reference outputs
  or tolerances were changed.
- `git diff --check` and Python syntax checks pass.

Evidence: `build-focused/xml-companion-contract-final-20261018/results.json`
and the individual case directories; `build-focused/Testing/Temporary/LastTest.log`.

## Real XML comparison

An initial NH3 periodic-stockholder HAR comparison converged through both the
existing library path and the new source-basis path. All 52 numerical fields
in the final refinement report were identical. Both gave R(F)=0.008335,
Rw(F)=0.007947, GoF=0.864031, 81 reflections, 13 fit parameters and scale
2.968414. This checks the same imported wavefunction through two basis routes;
it is not a comparison with a new theoretical calculation.

Evidence: `build-focused/xml-companion-real-nh3-20261017/comparison.json`.
The final executable was checked again with the strengthened real-data
script. Both fits converged and all 52 numerical fields were still identical.
All three reflection files were also **byte-identical** between the library
and companion paths:

- `xml_companion_real.archive.fcf`
- `xml_companion_real.archive.fco`
- `xml_companion_real.fcf6`

The final evidence is in
`build-focused/xml-companion-real-nh3-final-20261018/comparison.json`, with
input hashes and both runs retained in that directory. The source basis was
the explicit def2-TZVP block used for the original CRYSTAL calculation.

For this real XML the largest central-overlap difference is
`3.8317419801973829e-5`, within the existing `1e-4` acceptance threshold.
The machine-precision overlap result above belongs to the analytic fixture,
not to this external XML. The real test demonstrates equality of the source
and matching-library paths, not an exact zero overlap discrepancy.

The latest generic real-data checker also supports a repeatable, read-only
validation of the completed pair, which passed:

```sh
python3 scripts/check_crystal_xml_companion_real.py --check-existing \
  build-focused/xml-companion-real-nh3-final-20261018
```

It reports SHA-256 hashes for both copies of each identical output file. This
check does not rerun the fit, modify the supplied inputs or rewrite references.

## Limits and next action

This is opt-in and XML-only. It supports explicit general all-electron basis
blocks with S/SP/P/D/F/G shells and unit scale; it does not expand a named
internal basis, a full input file or an ECP basis. A central overlap check is
not a unique reconstruction of all radial functions. The supplied basis must
still come from the same calculation as the XML.

No MPI, macOS or non-GNU compiler run has been performed for this change.
No different private density or wavefunction interface was copied into this
port. The source and test scan found no such dependencies in the new files.

The original planned local checks passed. Additional investigation of Kang's
overlap reports is recorded below. No publication has been performed.

No pre-existing reference outputs have been changed and no commit or push
has been made.

## Follow-up: source-screened overlap entries

The exact companion basis does not eliminate differences caused by the
source calculation's integral screening. Failure diagnostics now distinguish
exact source zeros from retained overlaps, without changing either matrix
or the existing 1e-4 acceptance gate. New analytic tests check that accurate
source overlaps pass, deliberately zeroed overlaps above the gate fail, and
an incorrect basis still fails.

Public executable SHA-256 before the HA_fit guard described below:
`3f8fbff23db19af93d839c0c5a5185a69c931b885b685a986862d30ec403965b`.
The strengthened semantic checks pass 23/23; isolated short tests pass 73/73.
Full NH3 library/companion HAR comparisons remain identical in 52 numerical
report fields and all three reflection outputs. Against the previous build,
the reflection outputs differ only in the audit creation date. Existing
references and numerical tolerances have not been changed.

In the available L-alanine/def2-TZVP calculation, tightening source
TOLINTEG to `8 8 8 8 16` reduced the maximum central-overlap discrepancy
from 2.54732139778e-4 to 1.4731117455e-5. The source geometry was unchanged;
Tonto overlaps changed by less than 1.91e-13. Source SCF and XML export both
completed normally, and the existing check now accepts the regenerated
matrix. This confirms a screening remedy for this case, not a universal
bound or an automatic validation exception.

Those comparisons used central-block diagnostic fixtures. A full regenerated
XML HAR attempt exceeded the 900-second test-driver limit before refinement
output. This is incomplete validation, not a demonstrated numerical failure:
the original archived cycle took about 779 CPU seconds on its source machine.
A longer library trace identified an unintended molecular SCF: HA_fit
treated the missing finite-molecule density matrix as missing density even
though the XML periodic matrix was loaded separately. Added an oc-crystal23
guard to that initial-SCF fallback. No ANO/moment skipping, native GRED or
unrelated private feature was added. The redundant calculation was stopped
and both full HAR routes are being rerun. Independent full-XML density
comparisons pass on a coarse 5x5x5 grid with identical cube bytes; they
must not be presented as full HAR or electron-count validation. The HAR input
contains every measured amplitude/uncertainty retained in the archived FCF,
at printed precision, because the original HKL is unavailable locally. It
tests import-path equivalence, not exact historical refinement statistics.
The exact original Urea POB source calculation is also unavailable; workbook
integral checks must not be presented as reproduction of that export.

Reducing Tonto precision is not used as a remedy. The inspected XML overlap
tokens have 14 significant digits, while the workbook displays six decimal
places. The independently calculated Urea integral 0.000109530059260698
rounds to 0.000110 at six decimals, not zero. Coarse rounding can conceal a
mismatch but cannot recover screened source information.

After the HA_fit guard: rebuild passes and the isolated short suite passes
73/73. NH3 library/companion HAR converges with 52 identical report numbers;
all six reflection artifacts match the pre-guard run apart from audit date.
The real-data checker now explicitly exposes and rejects molecular-SCF
fallback while retaining the original numerical comparison tolerances.
Current public executable SHA-256:
`db21fea53740c5a55b56d012263194af73f83a9738aa08876f1eb12efa46e43d`.

## Completed follow-up, 21 October 2026

The full L-alanine library/companion HAR pair now passes: 283 final numeric
fields agree exactly and all three reflection files are byte-identical.
Public R(F)=0.018340, Rw(F)=0.015424, GoF=1.127228; 2373 reflections and
118 parameters. No molecular SCF fallback occurs. This remains one fixed-
density fit using archived-FCF amplitudes, not an outer HAR sequence or an
exact replay of the newly received Lolo_def2 first cycle.

The subsequently received Urea POB first-cycle archive is the exact source
of the supplied overlap spreadsheet: both complete printed 192-by-192
matrices match entry by entry. The unchanged gate rejects its four source-
zero discrepancies above 1e-4, with maximum 1.09530059261e-4. The same input
regenerated with TOLINTEG 8 8 8 8 16 and a matching restart converges and
exports complete XML. Geometry and Tonto overlaps are unchanged; maximum
discrepancy falls to 1.400286879e-6 and the four previously failing integrals
agree within 1.01e-15. This resolves that source-screening failure without
rounding, ignoring zeros or increasing the acceptance tolerance.

The complete regenerated XML then succeeds through both library and companion
HAR routes, using the received original CIF/HKL and settings. Both converge
in four LS steps with 817 selected reflections, 27 parameters, positive-
definite final ADPs and no molecular SCF. All 668 final numeric report fields
agree: R(F)=0.015892, Rw(F)=0.010255, GoF=2.458840, scale=3.150524. This is a
single fixed-density fit, not evidence that the final statistical model is
optimal or that a full outer HAR sequence has converged.

FCO outputs are byte-identical. The FCF files differ only in two exact
+180/-180 phase representations, at (0,0,-1) and (9,-9,0). The comparison
checker now recognizes exact whole-turn differences only in the explicitly
tagged `_refln_phase_calc` column. It reports byte identity separately and
retains hashes of the original files. No numerical tolerance is changed;
amplitudes, indices, row order, metadata and other fields must still agree.
Sixteen unit tests guard this comparison, including negative cases and
termination of the phase-column map at CIF data/save/stop boundaries. They are
registered as the short CTest `xml_companion_reflection_comparison`.

Received L-alanine/def2, Xylitol/def2 and Urea/def2 original XML failures also
consist entirely of above-threshold source zeros. Their exact source bases
do not remove that screening error. The two other newly failing Lolo_def2
transfers (NH3BF3 and hydrazine borane) remain incomplete; older passing
archives with different integral settings do not establish their resolution.
Do not present the example TOLINTEG setting as a universal guarantee.

No experimental input, large XML, private native-density reader or unrelated
development feature was added to this public worktree. No existing golden
output or numerical tolerance was loosened.

The separate exact Urea def2-TZVP source control has also completed. The same
geometry/basis with tighter screening reduces the maximum overlap discrepancy
from 2.26863475329e-4 to 4.0252355351e-6; Tonto's overlap matrix is bitwise
unchanged. All 16 originally failing entries agree within 4.90e-16 afterward,
and the source companion passes the unchanged gate in both builds. This is
an overlap/import test, not a completed def2 HAR refinement.

Final 21 October checks after comparator hardening: **74/74 public short
tests pass** (20.88 s), including all 16 comparison unit tests. The private
suite independently passes 95/95. All six retained full-XML NH3/L-alanine/
Urea-POB basis-route pairs pass read-only revalidation; final report numeric
differences are zero. No calculation output was rewritten for these checks.
The first public doctor test skipped because the non-login shell lacked the
already-installed mol2chemfig on PATH; rerunning with ~/.local/bin included
passed. Original logs are preserved. No scientific executable changed during
this final diagnostic/checker phase.

## GRED-free publication check, 21 October 2026

The production XML importer reads the requested XML and the explicitly
supplied general-basis companion. It does not read GRED or KRED to obtain
radial basis data, matrices or a fallback. An additional analytic contract
puts invalid GRED/KRED decoys in the run directory and requires unchanged
successful XML import. Public semantic checks now number 24; the private
variant adds one native-interface exclusion check (25).

New complete NH3 HAR calculations were traced with `strace -f -e trace=file`
in both builds. Their companion is byte-identical to the explicit basis
section in the original archived `.d12`; no native export was read to prepare
it. Both refinements converge, open XML and the companion, and make no
GRED/KRED file-access attempts. Structure-factor artifacts match the retained
results exactly except for their single creation-date field. File hashes
confirm that the original XML and basis were not modified. The legacy XML
reader's O_RDWR open mode is recorded, not misreported as O_RDONLY.

The original trace-harness failure is retained: it expected O_RDONLY and
byte identity across different creation dates. The separate successful
revalidation distinguishes those metadata checks from scientific results;
no general numerical or artifact-comparison tolerance was changed.

Final short-suite rerun including the native-file decoy case: **74/74 public**
(21.17 s) and **95/95 private** (52.48 s). No experimental data or traces
are included in the public repository. The Xylitol exact-source screening
control is separate continuing work; it does not gate these XML input checks
and is not claimed to have passed here.
