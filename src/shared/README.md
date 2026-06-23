# shared

Languages: English | [Korean](<README(ko).md>)

## What This Module Does

`shared` contains reusable building blocks for the Julia modules inside
`spectrum/`.

It is not a user-facing CLI. It provides common code for:

- crystal-cell and lattice representation,
- POSCAR/CONTCAR and `wannier90.win` lattice input,
- TOML parsing helpers with path resolution,
- Wannier `hr.dat` and `wsvec` data types and IO,
- Hermiticity completion and pair-symmetry checks for HR-like block dictionaries,
- `H(k)` construction and dense eigensystem solves,
- Wigner-Seitz vector generation and validation,
- Spglib-backed cell convention helpers used by structural supercells,
- projection group definitions from explicit indices, sidecar basis files, or
  `wannier90.win`,
- projected local atomic OAM operators for Wannier-basis observables,
- atomic angular momentum matrix conventions shared by OAM and SOC,
- spin angular momentum operators for spinful Wannier-basis observables,
- spinless-to-spinful HR and basis expansion helpers,
- and HR writing helpers used by supercell export.

## When To Use It

Read `shared` when you want to understand:

- how modules share a `CrystalCell`,
- how POSCAR/CONTCAR and `wannier90.win` lattice input is represented,
- how `hr.dat`, normalized hoppings, and `wsvec` tables are represented,
- how band, superham, and flake avoid duplicating TOML parsing utilities,
- how projection selectors are mapped to Wannier orbital indices,
- how projected atomic OAM matrices are assembled from orbital labels,
- how spin operators are assembled from complete `up/dn` basis pairs,
- how generated `wsvec` data is checked before it is used,
- how structural supercell generation standardizes POSCAR/CONTCAR inputs.

## Entry Points

There is no standalone CLI.

Load these helpers through a higher-level module such as:

```julia
include("spectrum/band/CLI.jl")
include("spectrum/superham/Core.jl")
```

Those module loaders include only the shared files they need.

## Key Files

- `CrystalCell.jl`: shared crystal-cell representation and lattice helpers.
- `PoscarIO.jl`: shared POSCAR reader built on `CrystalCell`.
- `LatticeIO.jl`: lattice readers for POSCAR/CONTCAR and `wannier90.win`
  `unit_cell_cart` blocks.
- `InputParsing.jl`: reusable TOML helper functions.
- `CellConventions.jl`: cell conventions, lattice transforms, and
  Spglib-backed standardized-cell helpers.
- `WannierTypes.jl`: shared types for HR and WS-vector data.
- `SpinLayout.jl`: shared `vasp544`/`qe` spin-layout parsing and index
  canonicalization.
- `WannierHrIO.jl`: shared `hr.dat` reading helpers.
- `Hermiticity.jl`: pair-Hermiticity error calculation and assertions for
  block dictionaries.
- `HrHermiticity.jl`: Hermitian partner completion and onsite block cleanup.
- `PairChecks.jl`: pair-dictionary validation helpers used by higher-level modules.
- `WannierKspace.jl`: plain and `wsvec`-aware Hamiltonian evaluation in `k`
  space.
- `WannierEigensystem.jl`: shared dense eigenvalue/eigenvector solves for Wannier Hamiltonians.
- `WannierWsvecIO.jl`: shared `wsvec` handling.
- `WannierWsvecGenerate.jl`: generated `wsvec` construction, validation, and
  writer utilities.
- `OrbitalProjection.jl`: projection group data model and eigenvector weights.
- `AtomicAngularMomentum.jl`: canonical p/d/t2g angular momentum matrices.
- `AtomicOam.jl`: local atomic `Lx`, `Ly`, `Lz`, and projected `L2` matrices.
- `AtomicSpin.jl`: spin `Sx`, `Sy`, `Sz`, and `S2` matrices for spinful bases.
- `SpinExpand.jl`: spinless HR/basis expansion and collinear up/down merging.
- `BasisLabelNormalize.jl`: shared orbital and spin label canonicalization.
- `Win90Basis.jl`: projection-basis reader and selector expansion for
  `wannier90.win`.
- `HrFormat.jl`: HR writing helpers shared by multiple modules.

## How It Works

`shared` does not run by itself.

Instead, other modules import pieces of it:

```text
band     -> LatticeIO, Hermiticity, PairChecks, WannierHrIO, WannierKspace, WannierEigensystem,
            SpinLayout, WannierWsvecIO, WannierWsvecGenerate, OrbitalProjection,
            AtomicOam, Win90Basis
soc      -> HrFormat, HrHermiticity, WannierHrIO, SpinLayout, SpinExpand,
            AtomicAngularMomentum, Win90Basis
superham -> HrFormat, CrystalCell, PoscarIO, SpinLayout, Hermiticity, PairChecks,
            WannierHrIO, WannierKspace, WannierEigensystem, WannierWsvecIO,
            WannierWsvecGenerate
supercell -> CrystalCell, CellConventions, PoscarIO, InputParsing
flake    -> shared helpers through superham/Core.jl
```

`InputParsing.namespaced_root(cfg, "band")` lets a module accept both
namespaced inputs such as `[band.run]` and standalone module-local inputs such
as `[run]`.

## Input And Output

There is no single input file for `shared`.

The lattice readers preserve the input cell basis:

```julia
include("spectrum/shared/CrystalCell.jl")
include("spectrum/shared/PoscarIO.jl")
include("spectrum/shared/LatticeIO.jl")

lattice = Main.LatticeIO.read_lattice("POSCAR")
```

`read_lattice` accepts POSCAR/CONTCAR-like files and `.win` files with a
`begin unit_cell_cart` block. TOML files are rejected for lattice input.

The path resolver used by module input parsers keeps relative paths relative to
the input TOML directory:

```julia
include("spectrum/shared/InputParsing.jl")

path = Main.InputParsing.resolve_path(dirname(abspath("input.toml")), "KPOINTS")
```

Projection helpers validate group labels, 1-based indices, duplicate coverage,
and whether a set of groups covers all Wannier orbitals.

`WannierWsvecGenerate.assert_wsvec_usable(hr, wsvec)` validates missing,
out-of-bounds, malformed, duplicate, and one-sided `wsvec` entries before a
workflow uses them.

## Dependencies

Internal dependencies:

- `CrystalCell.jl`
- `PoscarIO.jl`
- `LatticeIO.jl`
- `InputParsing.jl`
- `CellConventions.jl`
- `WannierTypes.jl`
- `SpinLayout.jl`
- `WannierHrIO.jl`
- `WannierKspace.jl`
- `WannierEigensystem.jl`
- `WannierWsvecIO.jl`
- `WannierWsvecGenerate.jl`
- `OrbitalProjection.jl`
- `AtomicOam.jl`
- `Win90Basis.jl`
- `HrFormat.jl`

External dependencies:

- Julia `LinearAlgebra`
- Julia `Printf`
- Julia `TOML`
- `Spglib`

## Minimal Code Snippet

Lattice-reading example:

```julia
include("spectrum/shared/CrystalCell.jl")
include("spectrum/shared/PoscarIO.jl")
include("spectrum/shared/LatticeIO.jl")

lattice = Main.LatticeIO.read_lattice("POSCAR")
println(lattice.source)
```

Projection example:

```julia
include("spectrum/band/CLI.jl")

group = Main.BandCLI.BandCore.OrbitalProjection.ProjectionGroup("A", [1, 2], "blue")
spec = Main.BandCLI.BandCore.OrbitalProjection.ProjectionSpec([group], 4)
```

`InputParsing.jl` is meant to be reused by module-specific `InputIO.jl` files,
not called as an application entry point.

## Run And Test

`shared` is tested indirectly through spectrum module checks.

Useful smoke checks:

```bash
julia --project=spectrum spectrum/main.jl band --help
julia --project=spectrum --startup-file=no spectrum/test/runtests.jl
```

## Common Pitfalls

- `spectrum` does not re-base structures internally; keep HR `R` vectors, reduced `KPOINTS`, and manual centers in the same input-cell basis.
- Fractional coordinates are meaningless unless you also know their reference cell.
- `run.wsvec` is explicit opt-in for band runs; merely placing a
  `wannier90_wsvec.dat` next to an input does not enable it.
- `wsvec` validation checks structural usability before numerical work starts,
  so incomplete tables fail early.
- `OrbitalProjection` uses 1-based Wannier indices, matching Julia and
  Wannier90-style orbital order.
- `AtomicOam` computes retained-basis projected local atomic OAM. It is not an
  orbital magnetization implementation.
- `AtomicSpin` requires complete spinful `up/dn` basis pairs; unpolarized basis
  metadata is rejected for SAM.

## Reading Order

1. `CrystalCell.jl`
2. `PoscarIO.jl`
3. `LatticeIO.jl`
4. `InputParsing.jl`
5. `CellConventions.jl`
6. `WannierTypes.jl`
7. `SpinLayout.jl`
8. `SpinExpand.jl`
9. `WannierHrIO.jl`
10. `HrHermiticity.jl`
11. `WannierKspace.jl`
12. `WannierEigensystem.jl`
13. `WannierWsvecIO.jl`
14. `WannierWsvecGenerate.jl`
15. `OrbitalProjection.jl`
16. `AtomicAngularMomentum.jl`
17. `AtomicOam.jl`
18. `Win90Basis.jl`
19. `HrFormat.jl`
