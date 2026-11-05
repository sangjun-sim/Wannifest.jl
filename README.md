# Wannifest.jl

## Introduction and structure

Wannifest.jl is a Julia package for Wannier tight-binding models. It works with Wannier90 Hamiltonians, POSCAR/CONTCAR structure inputs, Quantum ESPRESSO data, and TOML input files through a single command-line entry point.


| Command | Purpose |
| --- | --- |
| `band` | Band structures, DOS, projections, orbital and spin angular momentum |
| `contour` | Two-dimensional k-plane energy maps and contour/surface plots |
| `flux` | Directed complex flux terms, visualization, and diagnostics |
| `superham` | Supercell Hamiltonians from primitive-cell Wannier models |
| `soc` | SP/SOC Bloch overlaps and transformations between Wannier bases |


## Installation

Install Julia and make sure `julia` is on your `PATH`. The shell launcher also requires bash and `readlink -f`. From the repository root, install the dependencies listed in `Project.toml` and link the launcher:

```bash
julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate()'
mkdir -p "$HOME/.local/bin"
ln -s "$PWD/bin/wannifest" "$HOME/.local/bin/wannifest"
```

Skip the `ln` command if the symlink already points to this checkout. Add the
following to your shell configuration if `$HOME/.local/bin` is not on your `PATH`:

```bash
export PATH="$HOME/.local/bin:$PATH"
```


## Execution

From any working directory, show the available commands or run a TOML input:

```bash
wannifest --help
wannifest band input.toml
wannifest soc input.toml
```


Use `wannifest <command> --help` for command-specific options. The `--input input.toml` form is also supported. For TOML-backed commands, Wannifest can infer the command from the input path or top-level TOML namespace:


```bash
wannifest --input examples/band/rucl3/input.toml
wannifest examples/contour/graphene/input.toml
```

Command-line input paths are relative to the current working directory. Paths inside TOML inputs are generally relative to the input file; output locations are controlled by the command's input settings and command-line options.

## Regenerating the Fe SOC example

The large QE output directories are excluded from Git. The example retains
`scf.in`, `nscf.in`, the pw2wannier90 inputs, and the three Wannier90 `.win` files.
QE and Wannier90 are needed to regenerate the data; they are not needed for the
automated tests. The QE inputs use `pseudo_dir = './pseudo'`. Place the
scalar-relativistic ONCV `Fe.upf` in `sp/pseudo/` and the fully relativistic ONCV
`Fe.upf` in `soc/pseudo/`; these pseudopotentials are not included.

Run these commands in `examples/soc/fe/sp/`:

```bash
pw.x < scf.in > scf.out
pw.x < nscf.in > nscf.out
wannier90.x -pp fe_up
wannier90.x -pp fe_dn
pw2wannier90.x < pw2wan_up.in > pw2wan_up.out
pw2wannier90.x < pw2wan_dn.in > pw2wan_dn.out
wannier90.x fe_up
wannier90.x fe_dn
```

Then run these commands in `examples/soc/fe/soc/`:

```bash
pw.x < scf.in > scf.out
pw.x < nscf.in > nscf.out
wannier90.x -pp fe
pw2wannier90.x < pw2wan.in > pw2wan.out
wannier90.x fe
```

Keep the non-HDF5 wavefunction files and `data-file-schema.xml` in each channel's
`outputs/fe.save/`, as configured in `input.toml`. Recreate all three Wannier90
channels after rerunning QE: the supplied `.eig`, `_u.mat`, and `_u_dis.mat`
files must correspond to the same NSCF wavefunctions. An SCF input alone does
not replace those outputs. From the repository root, run:

```bash
wannifest examples/soc/fe/input.toml
```

## Tests

From the repository root, run the test suite with the dependencies in
`Project.toml` and Julia's `Test` standard library:

```bash
GKSwstype=100 JULIA_NUM_THREADS=2 julia --project=. -e 'using Pkg; Pkg.test()'
```

GitHub Actions runs the same suite on pushes and pull requests with Julia 1.12.
The tests check command inference, example input parsing, analytic band and
contour energies, degeneracy and wsvec weights, flux continuity, supercell band
folding, orbital/spin angular momentum, and SOC k-point matching, Fourier phases,
and binary output validation. Synthetic inputs and calculation outputs are
created in system temporary directories and removed after the tests. No display
or external QE calculation data are required; full QE/SOC calculations and GUI
behavior are outside this suite.
