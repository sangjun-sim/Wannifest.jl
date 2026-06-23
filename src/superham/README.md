# superham

Languages: English | [Korean](<README(ko).md>)

## What This Module Does

`superham` builds a supercell Hamiltonian from a primitive-cell `hr.dat`.
The structure lattice is used exactly as provided; `hr.dat` `R` vectors and any
manual centers must already be expressed in that same input-cell basis.

It can also attach geometry information and Wigner-Seitz corrections from:

- `wsvec` tables,
- Wannier centers from `centres.xyz`,
- user-provided manual centers,
- or atom-based surrogate centers.

The CLI can optionally check a folded spectrum at a requested `k` point and
write a normalized supercell `hr.dat`.

## When To Use It

Use `superham` when you want to:

- expand a primitive Hamiltonian into a larger supercell Hamiltonian,
- keep or attach center information for geometry-aware workflows,
- test how a tight-binding model changes under supercell folding,
- write a new `hr.dat` for a supercell.

## Entry Points

Main entry:

- `../main.jl`

Module CLI:

- `CLI.jl`

Core implementation:

- `Service.jl`

## Key Files

- `CLI.jl`: command-line interface.
- `InputIO.jl`: reads files, geometry, and supercell options from TOML.
- `Service.jl`: main orchestration.
- `HrIO.jl`: reads the input `hr.dat`.
- `CenterIO.jl`: loads or builds center tables.
- `SupercellHam.jl`: constructs the expanded Hamiltonian.
- `Kspace.jl`: evaluates eigenvalues at `k`.
- `Validate.jl`: Hermiticity and folding checks.
- `WsvecIO.jl`: optional Wigner-Seitz correction support and generated
  `wsvec` entry point.
- `../shared/HrFormat.jl`: writes normalized supercell `hr.dat` output.
- `../test/runtests.jl`: shared spectrum tests that cover superham helpers.

## How It Works

High-level flow:

```text
CLI -> InputIO -> Service -> PoscarIO.read_poscar
                         -> HrIO.read_hr
                         -> optional WsvecIO.attach_wsvec
                         -> optional CenterIO.centers_from_config
                         -> SupercellHam.build_supercell_model
                         -> Validate hermiticity / size / folded spectrum
                         -> optional normalized hr.dat writer
```

The module supports several geometry modes:

- `none`: no orbital-center table is attached.
- `wsvec`: signals a wsvec-only geometry policy; actual `wsvec` attachment is
  controlled by the non-empty `files.wsvec` path.
- `manual_centers`: use the `[[superham.geometry.orbitals]]` entries as
  fractional Wannier centers.
- `atomic_assumption`: use `atom_index` in each orbital entry to copy atom
  positions from the structure file as approximate centers.

The geometry mode controls whether orbital-center information is omitted,
provided manually, approximated from atoms, or paired with explicit `wsvec`
data.

## Input And Output

Minimal input example:

```toml
[superham.files]
hr = "wannier90_hr.dat"
structure = "POSCAR"
win = ""
wsvec = ""
centres = ""
output_hr = "wannier90_sc_hr.dat"

[superham.geometry]
mode = "none"
strict = false

[superham.supercell]
matrix = [2, 0, 0, 0, 2, 0, 0, 0, 1]
```

The supercell matrix may also be written as a nested 3x3 array:

```toml
[superham.supercell]
matrix = [
  [2, 0, 0],
  [0, 2, 0],
  [0, 0, 1],
]
```

Manual-centers example:

```toml
[superham.geometry]
mode = "manual_centers"
manual_num_wann = 2

[[superham.geometry.orbitals]]
label = "A:pz"
center_frac = [0.0, 0.0, 0.0]

[[superham.geometry.orbitals]]
label = "B:pz"
center_frac = [0.3333333333, 0.3333333333, 0.0]
```

Atomic-assumption example:

```toml
[superham.geometry]
mode = "atomic_assumption"
manual_num_wann = 2

[[superham.geometry.orbitals]]
label = "A:pz"
center_frac = [0.0, 0.0, 0.0]
atom_index = 1

[[superham.geometry.orbitals]]
label = "B:pz"
center_frac = [0.0, 0.0, 0.0]
atom_index = 2
```

Run example:

```bash
julia main.jl superham --input input.superham.toml
julia main.jl superham --input input.superham.toml --kpoint 0 0 0 --output-hr wannier90_sc_hr.dat
```

Outputs:

- supercell `hr.dat`,
- stdout summary with geometry source, primitive/supercell orbital counts,
  folded-spectrum validation, and eigenvalues at the selected `k` point.

`output_hr` in `[superham.files]` is used unless the CLI `--output-hr` override
is provided. If neither is set, the service still runs validations and prints a
summary without writing an HR file.

`files.win` is accepted in the configuration object, but the current service
builds its lattice and atom positions from `files.structure`.

## Dependencies

Internal dependencies:

- `../shared/CrystalCell.jl`
- `../shared/HrFormat.jl`
- `../shared/WannierHrIO.jl`
- `../shared/WannierKspace.jl`
- `../shared/WannierWsvecIO.jl`
- `../shared/WannierWsvecGenerate.jl`
- `../shared/PoscarIO.jl`
- `InputIO.jl`
- `CenterIO.jl`
- `SupercellHam.jl`
- `Validate.jl`

External dependencies:

- Julia `TOML`
- Julia `LinearAlgebra`

## Minimal Code Snippet

CLI usage:

```bash
julia main.jl superham --input input.superham.toml
```

Service-layer usage:

```julia
include("spectrum/superham/CLI.jl")

cfg = Main.SuperhamCLI.InputIO.read_input("input.superham.toml")
result = Main.SuperhamCLI.Service.run(cfg)
println(result.output_hr)
```

## Run And Test

Run the module:

```bash
julia main.jl superham --input input.superham.toml
```

Run tests:

```bash
julia --project=spectrum --startup-file=no spectrum/test/runtests.jl
```

## Common Pitfalls

- Structure, HR `R` vectors, and manual center coordinates must use the same input-cell basis.
- `geometry.mode = "manual_centers"` and `geometry.mode = "atomic_assumption"` require `geometry.orbitals`.
- `manual_num_wann` must match `hr.dat` basis size.
- `atomic_assumption` is approximate unless the Wannier centers are actually
  atom-centered; the code warns when it uses this mode.
- `centres.xyz` takes precedence over manual or atomic center generation when
  the file path is set and exists.
- `files.wsvec = ""` disables `wsvec`; a non-empty path must pass coverage
  validation before it is attached.
- `geometry.strict = true` makes geometry-sensitive supercell construction more
  conservative about missing center information.

## Reading Order

1. `CLI.jl`
2. `InputIO.jl`
3. `Service.jl`
4. `HrIO.jl`
5. `CenterIO.jl`
6. `SupercellHam.jl`
7. `Validate.jl`
8. `../test/runtests.jl`
