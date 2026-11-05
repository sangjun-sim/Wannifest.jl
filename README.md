# Wannifest.jl

## Introduction and structure

Wannifest.jl is a Julia package for Wannier tight-binding models. It works with Wannier90 Hamiltonians, POSCAR/CONTCAR structure inputs, Quantum ESPRESSO data, and TOML input files through a single command-line entry point.


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
