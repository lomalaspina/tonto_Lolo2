# Reading the source basis with CRYSTAL23 XML

The CRYSTAL23 XML file contains the density matrix and, when exported, overlap
and Fock matrices. It does not contain the primitive Gaussian exponents and
contraction coefficients. Tonto can read those from a separate, explicit
CRYSTAL basis block from the **same calculation**.

This is an optional XML-only input. Without it, the existing `basis_name` and
`basis_directory` workflow is unchanged. The reader does not derive a basis
from the overlap matrix or adjust the density to force agreement.

The XML route does not open GRED or KRED, including to obtain the basis.
Supply the explicit basis block used in the original CRYSTAL input (or the
same full-precision source definition), not a native wavefunction export.
There is no automatic native-file fallback or companion-generation step.
If the original input selects an internal basis by name, obtain its exact
explicit definition or retain the matching Tonto library workflow; the name
alone and the inspected XML cannot supply missing radial coefficients.

## Input

```text
{
   name= crystal_example
   CIF= { file_name= crystal_example.cif }

   ! Still needed for the Thakkar proatoms used in Hirshfeld partitioning.
   basis_directory= /path/to/tonto/basis_sets

   c23_XML_file_name= GenerateXML.XML
   c23_xml_basis_file= exact.crystal.basis
   process_cif_and_c23_xml

   ! Continue with the usual crystal and refinement settings.
}
```

`exact.crystal.basis` must contain only the explicit general-basis section,
including its final `99 0` line. Do not give a complete `.d12` input or a
printed output listing. Preserve the exponents and coefficients at their
original precision. A basis-set name is not enough.

The imported Gaussian basis does not need a Tonto `basis_name`. Hirshfeld
partitioning still needs Tonto's spherical proatom data (normally Thakkar),
so keep a valid `basis_directory` when continuing to a HAR calculation.

For example, the following shows the syntax of a hydrogen basis block:

```text
1 1
0 0 3 1.0 1.0
3.425250914 0.1543289673
0.6239137298 0.5353281423
0.1688554040 0.4446345422
99 0
```

This syntax example is not a recommended basis for a periodic calculation.
Use the exact basis used to generate the accompanying XML.

## Supported shells and atom assignment

The record definitions follow section 2.2.1 of the
[CRYSTAL23 User's Manual](https://www.crystal.unito.it/include/manuals/crystal23.pdf).
The companion reader supports the following explicit subset of that format.

- General all-electron shells, `ITYB=0`, are supported.
- `LAT=0,1,2,3,4,5` means S, SP, P, D, F and G, respectively.
- SP shells are split into their S and P components without changing their
  exponents or coefficients.
- The scale field must be `SCAL=1`. Internal/standard basis codes,
  pseudopotentials and other scale values are rejected.
- Blank lines and `!` or `#` comments are allowed. Extra fields, invalid or
  non-finite numbers and zero contractions are rejected.

One unambiguous block per element is assigned automatically to all XML atoms
of that element. Definition blocks need not follow the order of the atoms.
If different atoms of the same element use different bases, use distinct
CRYSTAL atomic identifiers and append an explicit atom map after `99 0`:

```text
ATOM_MAP 2
1 1
2 101
```

This maps XML atom 1 to the block with identifier 1 and XML atom 2 to block
101. Both identifiers denote hydrogen. Supply one row for **every atom in
the original XML**, before Tonto sorts atoms. Each identifier must refer to
a supplied block of the correct element. Ambiguous assignments are rejected.

## Checks and limitations

The import keeps an individual basis for each XML atom, checks the number
and component order of its spherical AOs, and applies Tonto's Gaussian
normalization conversion. XML still supplies the matrices; their existing
atom/AO permutation and lattice-block orientation are unchanged.

The XML must contain its central-cell overlap matrix. Tonto constructs that
matrix from the imported basis and requires agreement within `1e-4` in the
largest matrix element difference. A wrong basis can therefore be rejected
even when its AO count is correct. The check runs before the density is
used for partitioning. A missing overlap matrix is an error for this option.

Central overlap agreement is a necessary check, not proof that two basis
definitions are identical at every point in space. The user must keep the
basis and XML from the same source calculation. This option does not repair
incomplete or inconsistent XML, infer missing radial functions, or establish
that a chosen molecular basis is well conditioned in a periodic system.

### Avoiding repeated overlap integrals

After a successful check, Tonto retains an in-memory snapshot for that XML
import. A later check reuses the successful result only when the actual basis
exponents/contractions, atom positions, shell/AO layout, ordering settings and
imported central overlap values are exactly unchanged. Comparison uses no
rounding, hash or numerical tolerance. Dimension/mapping and finite-value
safeguards remain active; only repeated Gaussian-overlap evaluation is skipped.

Every new XML read starts fresh, even if its filename and contents are the
same. Rebuilding the shell-pair data also invalidates the snapshot. Geometry,
basis or source-overlap changes require a fresh successful comparison. The
cache belongs to the import object and is not saved between jobs; passing one
HAR cycle does not exempt subsequent XML imports from validation.

There is no user bypass switch and the overlap tolerance remains `1e-4`.
This optimization does not change CRYSTAL's `TOLINTEG`, any imported matrix,
the reconstructed density or the refinement model. It requires one additional
real central-overlap matrix in memory (approximately `8*nAO*nAO` bytes).

The short test `xml_overlap_cache_integration` loads the analytic XML fixture
through the production reader, then checks cache hits and invalidation after
rereads and changes to geometry, basis, overlap or ordering. Invalid dimensions,
atom/AO mappings and NaN/Inf values must be rejected even after a successful
cached check. It also verifies that a hit does not modify the imported overlap
or density matrices. Build `run_crystal_xml_overlap_cache`, then run:

```sh
ctest --test-dir <build-directory> -R '^xml_overlap_cache_integration$' --output-on-failure
```

CRYSTAL's exported overlap can contain screened-out entries, whereas Tonto
evaluates the corresponding Gaussian overlaps explicitly. Consequently, the
current maximum-element check can reject a correct source basis when a small
nonzero Tonto overlap corresponds to an exported zero. Importing the exact
companion basis does not remove this screening discrepancy. Inspect the
source integral thresholds and the differing entries before interpreting
the failure as a basis mismatch. Do not simply ignore every zero entry or
increase the acceptance tolerance: either can hide a genuine basis/order
error. On failure, Tonto reports separate discrepancies and counts for exact
source zeros and retained overlaps. This is diagnostic information, not an
automatic exception to the acceptance rule.

For screening-related discrepancies, regenerate the CRYSTAL calculation and
its XML with tighter integral settings, keeping the basis, geometry and
k-point mesh consistent. For example, the SCF block can include:

```text
TOLINTEG
8 8 8 8 16
```

These are integral-screening controls, not decimal-output controls. Printing
more digits cannot recover a screened integral. Check convergence and the
new overlap comparison; this example setting is not a universal accuracy
guarantee. Keep the original results and do not substitute only one matrix
from an unrelated calculation.

In a controlled Urea/pob-TZVP-rev2 example, the exact source basis already
matched Tonto, but four exported zero entries differed by more than `1e-4`.
Regenerating the same geometry/basis with the example settings reduced the
largest discrepancy from `1.0953e-4` to `1.4003e-6`. Tonto's overlap matrix
was bitwise unchanged, and the four formerly missing entries then agreed
within `1.01e-15`. The regenerated complete XML also passed a fixed-density
HAR comparison between the source-companion and matching-library routes.
This is evidence for a source-screening remedy in that case, not for
reducing Tonto precision or automatically accepting arbitrary screened zeros.

An independent Urea/def2-TZVP control showed the same behavior: the maximum
difference fell from `2.2686e-4` to `4.0252e-6`, with Tonto's overlap matrix
unchanged and all 16 formerly failing entries agreeing within `4.90e-16`.
The regenerated exact-source companion imports pass in both tested builds.
This latter check validates overlaps/import, not a full def2 HAR sequence.

Do not reduce Tonto's computational precision to match a rounded matrix
listing. A listing or spreadsheet may show fewer digits than the XML actually
stores. Rounding to six decimal places changes an entry by at most `5e-7`;
it cannot explain an exported zero whose independently calculated overlap is
about `1e-4`. Rounding both matrices more coarsely can conceal a mismatch
without restoring the missing source information. These considerations also
apply when the basis was imported using `c23_xml_basis_file`.

See the [short-test description](../tests/short/crystal_xml_companion_basis/README.md)
for the analytic accepted/rejected-input checks and their limits.
