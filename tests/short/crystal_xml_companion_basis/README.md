# XML companion-basis import contract

This is an **analytic interface fixture**, not a CRYSTAL calculation or a
refinement benchmark. It is deliberately small enough for the short suite.
The two inequivalent H sites have normalized primitive s Gaussians with
exponents 0.7 and 1.3 bohr^-2, separated by 2 bohr in a 10-bohr P1 cell.
A He site uses exponent 1.1 bohr^-2 and is 2 bohr from H1 and sqrt(8)
bohr from H2. The source XML order is H1/He/H2, so the importer must also
reorder atom/AO blocks while retaining the atom-specific basis assignment.
The exact pairwise central overlap is

`S12 = (4*a*b/(a+b)^2)^(3/4) * exp[-a*b*R^2/(a+b)]`

The H1-H2 overlap is `0.15096111232544537`, H1-He is
`0.17393146951240948`, and He-H2 is `0.0084643405576893919`.
Diagonal overlaps are one. The illustrative source density is diag(1,2,1)
(electron count `Tr(PS)=4`); a symmetric illustrative Fock matrix is included.
Only R=0 is supplied.
These values are independently analytic and are not fitted to Tonto output.

`source.basis` contains explicit CRYSTAL general-shell data and a source-atom
map. NAT identifiers 1 and 101 denote different bases on the same element.
The input intentionally sets no `basis_name` or `basis_directory`: importing
the companion must provide all required basis information. Definitions are
also tested in reversed order to ensure that element/block order does not
replace the source-atom map.

The fixture and its companion are analytic and require no GRED or KRED file.
A separate control adds deliberately invalid native-export decoys to the
job directory; the XML import must still succeed unchanged. Those unrelated
files must not supply a hidden basis or matrix fallback.

The semantic checker also swaps the assignments of the two H basis sets:
this must fail their unequal overlaps with He even though both H atoms have
one AO. It runs the valid import plus deliberately invalid basis,
AO-count, overlap, map, exponent, scale and AO-label controls. A segmentation
fault is never an accepted negative-test result. Same-size wrong exponents
must fail the overlap check before any imported density is used.
Nonfinite NaN/Infinity entries in each of the density, overlap and Fock
matrices must also be rejected; a successful overlap test cannot establish
that unrelated matrix entries are finite.
A same-process second import without the optional Fock matrix must report
that Fock data are absent, rather than reusing the previous file's matrix.
Similarly, a second import without the required central overlap must fail.

An additional analytic pair has overlap 0.00025. Its accurate source value
passes, while replacing it with zero still fails the unchanged tolerance.
The failure must report that the discrepancy is at source-zero entries and
explain how to investigate source screening. A wrong basis with the accurate
small-overlap source must also fail. This tests the diagnostic and protects
against silently accepting every zero as a screened integral.

Run with:

```sh
python3 scripts/check_crystal_xml_companion_basis.py \
  --program build/tonto \
  --fixture tests/short/crystal_xml_companion_basis
```

Each run retains its own inputs, stdout and driver log and reports the path.
No reference is automatically blessed. This test alone does not establish
real CRYSTAL density correctness. The separately recorded source-basis NH3
comparison exercises full HAR with actual source exponents/contractions.
The existing canonical ammonium-borane test continues to cover the legacy
library-based XML workflow; it is not itself a companion-basis benchmark.
