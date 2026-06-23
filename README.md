# Wannifest.jl

Wannifest.jl is a Julia toolkit for working with Wannier tight-binding models
and crystal-structure inputs. It provides a single command-line entry point for
band analysis, k-plane energy maps, flux-term construction, structural
supercells, and supercell Hamiltonians.

The package is designed for workflows built around Wannier90-style `hr.dat`
files, `KPOINTS`, POSCAR/CONTCAR structures, `wannier90.win` projection data,
and TOML input files. It emphasizes explicit input schemas, reproducible
generated data, and command-level tools that can be run independently.

## Capabilities

- `band`: band structure, DOS, projection, OAM, and SAM workflows.
- `contour`: two-dimensional k-plane energy surfaces and contour/surface plots.
- `flux`: directed complex flux terms for Wannier Hamiltonians, with optional
  3D HTML visualization and diagnostics.
- `supercell`: POSCAR/CONTCAR structural supercell generation with symmetry
  summaries and validation output.
- `superham`: primitive-to-supercell Hamiltonian construction with optional
  orbital-center and Wigner-Seitz vector handling.

## Installation

Instantiate the Julia environment from the repository root:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Usage

Show the command list:

```bash
julia --project=. main.jl --help
```

Run a band calculation without plotting:

```bash
julia --project=. main.jl band --input examples/band/rucl3/input.toml --no-plot
```

Sample a k-plane energy surface:

```bash
julia --project=. main.jl contour --input examples/contour/graphene/input.toml --no-plot
```

Build a structural supercell from `POSCAR` in the current directory:

```bash
julia --project=. main.jl supercell 2 0 0, 0 2 0, 0 0 1
```

Add directed flux terms and write only the Hamiltonian output:

```bash
julia --project=. main.jl flux --input examples/flux/graphene/input.toml --no-html --no-diagnostic
```

Build a supercell Hamiltonian from a TOML input:

```bash
julia --project=. main.jl superham --input path/to/input.superham.toml
```

For TOML-backed commands, the command can also be inferred from the input path
or top-level TOML namespace:

```bash
julia --project=. main.jl --input examples/band/rucl3/input.toml
julia --project=. main.jl examples/contour/graphene/input.toml
```

## Inputs And Outputs

Most commands use namespaced TOML tables such as `[band.run]`,
`[contour.plane]`, `[flux.terms]`, or `[superham.supercell]`. Paths are resolved
relative to the input file unless a command-line override is provided.

Generated data and plot files are written under command-specific `outputs/`
directories by default. Example outputs are ignored by git so examples can be
rerun without polluting the working tree.

## Repository Layout

- `main.jl`: executable command-line entry point.
- `src/Wannifest.jl`: command dispatch and TOML command inference.
- `src/shared/`: shared Wannier, lattice, spin-layout, projection, and parsing
  utilities.
- `src/band/`, `src/contour/`, `src/flux/`, `src/supercell/`,
  `src/superham/`: command implementations.
- `examples/`: runnable example inputs and reference data.
- `test/`: integration and module tests.

## Tests

```bash
julia --project=. --startup-file=no test/runtests.jl
```
