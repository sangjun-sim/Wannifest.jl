# supercell

Languages: English | [Korean](<README(ko).md>)

## What This Module Does

`supercell` builds structural supercells from POSCAR/CONTCAR structure input.

It applies an integer supercell matrix directly to the input basis, writes a
new POSCAR, prints Spglib symmetry summaries, and validates symmetry mappings.

## When To Use It

Use `supercell` when you want to:

- expand a structure into a larger POSCAR,
- check the symmetry summary before and after supercell construction,
- write mismatch diagnostics for atom-symmetry mapping checks.

## Entry Points

Main entry:

- `../main.jl` (`spectrum/main.jl`)

Module CLI:

- `CLI.jl`

Core implementation:

- `Service.jl`

## Key Files

- `CLI.jl`: parses arguments and prints symmetry/check summaries.
- `InputIO.jl`: builds `RunConfig` from CLI matrix/path arguments and fixed
  command defaults.
- `Service.jl`: loads input structures, builds the supercell in the original
  input basis, writes POSCAR output, and runs optional validation.
- `Core.jl`: loads the supercell submodules.
- `Symmetry.jl`: Spglib-backed dataset summary helpers.
- `Transform.jl`: integer-matrix supercell construction.
- `Validate.jl`: optional symmetry-mapping validation and mismatch output.

## How It Works

High-level flow:

```text
CLI -> InputIO -> Service -> load POSCAR/CONTCAR structure
                         -> optional Spglib dataset summary
                         -> Transform.build_supercell
                         -> write POSCAR
                         -> Validate checks
```

The structural basis is always the original input POSCAR/CONTCAR basis.
The CLI fixes `use_symmetry=true`, `validate=true`, `symprec=1.0e-5`,
`angle_tolerance=-1.0`, and `digits=12`.

## Input And Output

Minimal CLI examples:

```bash
julia main.jl supercell 2 2 1
julia main.jl supercell 2 0 0, 0 2 0, 0 0 1
```

The command reads `POSCAR`, writes `POSCAR.supercell`, and
writes the mismatch table next to the output POSCAR as `atom_symm_mismatch.dat`.

Outputs:

- supercell POSCAR,
- stdout symmetry/check summary,
- optional atom-symmetry mismatch table.

## Dependencies

Internal dependencies:

- `../shared/CrystalCell.jl`
- `../shared/CellConventions.jl`
- `../shared/PoscarIO.jl`
- `Core.jl`
- `Symmetry.jl`
- `Transform.jl`
- `Validate.jl`

External dependencies:

- `Spglib`
- Julia `LinearAlgebra`

## Minimal Code Snippet

CLI usage:

```bash
julia main.jl supercell 2 2 1
```

Service-layer usage:

```julia
include("spectrum/supercell/CLI.jl")

matrix = Main.SupercellCLI.InputIO.parse_matrix_args(["2", "2", "1"])
cfg = Main.SupercellCLI.InputIO.default_config(matrix)
result = Main.SupercellCLI.Service.run(cfg)
println(result.output_path)
```

## Run And Test

Run the module:

```bash
julia main.jl supercell 2 2 1
```

Smoke-check the CLI:

```bash
julia --project=spectrum spectrum/main.jl supercell --help
```

## Common Pitfalls

- `supercell` does not accept `--input`, `--output`, or other extra arguments.
- The command always reads `POSCAR`; if that file does not exist, it errors.
- Use commas, not semicolons, between full-matrix rows.
- Three matrix entries are interpreted as diagonal values. A full matrix must be
  three comma-separated rows with three integer entries each.
- Validation is always enabled from the CLI and writes a mismatch table.

## Reading Order

1. `CLI.jl`
2. `InputIO.jl`
3. `Service.jl`
4. `Symmetry.jl`
5. `Transform.jl`
6. `Validate.jl`
