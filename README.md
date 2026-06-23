# Wannifest.jl

Wannifest is a copy-only Julia package containing selected spectrum workflows:

- `band`: band structure, DOS, projection, OAM, and SAM workflows.
- `contour`: 2D k-plane energy surface data and plots.
- `flux`: directed flux terms for Wannier Hamiltonians.
- `supercell`: structural POSCAR/CONTCAR supercell builder.
- `superham`: supercell Hamiltonian builder.

The source repository remains intact. This package intentionally duplicates the
shared helper layer needed by the copied commands.

## Usage

```bash
julia --project=/Users/sangjun/Wannifest.jl /Users/sangjun/Wannifest.jl/main.jl --help
julia --project=/Users/sangjun/Wannifest.jl /Users/sangjun/Wannifest.jl/main.jl band --input /Users/sangjun/Wannifest.jl/examples/flake/graphene/input.band.toml --no-plot
julia --project=/Users/sangjun/Wannifest.jl /Users/sangjun/Wannifest.jl/main.jl contour --input /Users/sangjun/Wannifest.jl/examples/contour/graphene/input.toml --no-plot
julia --project=/Users/sangjun/Wannifest.jl /Users/sangjun/Wannifest.jl/main.jl superham --input /Users/sangjun/Wannifest.jl/examples/flake/graphene/input.superham.toml
julia --project=/Users/sangjun/Wannifest.jl /Users/sangjun/Wannifest.jl/main.jl supercell 2 0 0, 0 2 0, 0 0 1
julia --project=/Users/sangjun/Wannifest.jl /Users/sangjun/Wannifest.jl/main.jl flux --input /Users/sangjun/Wannifest.jl/examples/flux/graphene/input.toml --no-html --no-diagnostic
```

## Not Included

This copy does not include flake, SOC, or utility commands.

## Tests

```bash
julia --project=/Users/sangjun/Wannifest.jl --startup-file=no /Users/sangjun/Wannifest.jl/test/runtests.jl
```
