# band

Languages: English | [Korean](<README(ko).md>)

## What This Module Does

`band` computes:

- band structures,
- density of states (DOS),
- orbital-projected band weights and projected DOS (PDOS),
- projected local atomic OAM along a band path,
- spin angular momentum expectation along a band path,
- spin-colored plots and spin-resolved DOS when an explicit spin layout is
  requested,
- or any supported combination of those outputs,

from a Wannier-style `hr.dat`.

It supports explicit opt-in `wsvec` corrections, POSCAR/CONTCAR or
`wannier90.win` lattice input, lazy plotting through `Plots`, and strict TOML
validation.

## When To Use It

Use `band` when you want to:

- plot a band path from a `KPOINTS` file,
- compute a DOS curve on a uniform mesh,
- inspect orbital character along a band path or in DOS,
- compare bands from different Hamiltonians,
- produce data files for external plotting.

## Entry Points

Main entry:

- `../main.jl`

Module CLI:

- `CLI.jl`

Core implementation:

- `Service.jl`

## Key Files

- `CLI.jl`: parses `--input` and `--no-plot`.
- `InputIO.jl`: reads TOML and validates run, output, DOS, plotting, spin,
  projection, and OAM options.
- `Service.jl`: main orchestration for band and DOS runs.
- `../shared/LatticeIO.jl`: reads POSCAR/CONTCAR or `wannier90.win` lattice input.
- `KPath.jl`: parses `KPOINTS` and generates interpolated paths.
- `Execution.jl`: temporarily sets BLAS thread counts for scoped execution.
- `Bands.jl`: band-structure post-processing.
- `Dos.jl`: DOS mesh generation and broadening.
- `projection/Model.jl`: projection config and result types.
- `projection/Input.jl`: `[band.projection]` parser and default output paths.
- `projection/InputGroups.jl`: compact/table projection group parser.
- `PlotInput.jl`: `[band.plot]` parser and explicit plot-target validation.
- `projection/Core.jl`: builds projection specs and projection outputs.
- `observables/Oam.jl`: builds projected atomic OAM contexts from projection
  basis metadata.
- `observables/Sam.jl`: builds spin angular momentum contexts from spinful basis
  metadata.
- `observables/Output.jl`: writes long-format OAM and SAM data.
- `observables/Plot.jl`: plots OAM- and SAM-colored band structures.
- `../shared/SpinLayout.jl`: shared explicit spin-layout indices, spin
  coloring, and spin summary text.
- `Output.jl`: data writers for bands, DOS, projection weights, and PDOS.
- `PlotService.jl`: lazy loader and dispatcher for plotting.
- `PlotBands.jl`: plotting backend through `Plots`, loaded only when plotting
  is requested.
- `../shared/OrbitalProjection.jl`: group definitions and eigenvector weight calculation.
- `../shared/AtomicOam.jl`: local atomic `Lx`, `Ly`, `Lz`, and `L2` operators.
- `../shared/Win90Basis.jl`: reads `wannier90.win` projections for selector-based groups.
- `../shared/WannierEigensystem.jl`: diagonalizes `H(k)`.

## How It Works

High-level flow:

```text
CLI -> InputIO -> Service -> layout-aware WannierHrIO.read_hr
                         -> optional layout-aware WannierWsvecIO.read_wsvec + validation
                         -> optional LatticeIO.read_lattice
                         -> KPath / Dos mesh generation
                         -> WannierEigensystem.solve_kpoint
                         -> explicit spin layout split when enabled
                         -> optional projection weights / PDOS
                         -> optional projected atomic OAM
                         -> Output data writers
                         -> optional PlotService plots
```

There are three common compute paths:

- band path: parse `KPOINTS`, interpolate a path, solve eigenvalues, write band data,
- DOS path: build a regular mesh, solve many `k` points, broaden eigenvalues, write DOS data.
- projection path: use sorted eigenvectors from the same diagonalization, sum `abs2` over group indices, and write separate projection outputs.
- OAM path: use sorted eigenvectors from the same diagonalization and write
  projected local atomic OAM as separate long-format data.

## Input And Output

Minimal band example:

```toml
[band.run]
mode = "bands"
# Defaults when omitted:
# hr = "wannier90_hr.dat"
# wsvec omitted or wsvec = "" disables wsvec corrections.
# Set wsvec = "wannier90_wsvec.dat" explicitly to enable them.
# kpoints = "KPOINTS"
# structure = "POSCAR"
# hermiticity_tol = 1.0e-8
# verbose = true

[band.energy]
shift = 0.0
```

Minimal DOS example:

```toml
[band.run]
mode = "dos"
hr = "wannier90_hr.dat"
structure = "POSCAR"

[band.dos]
mesh = "50x50x50"
shift = [0.0, 0.0, 0.0]
sigma = 0.05
npts = 2001
emin = -10
emax = 10

[band.energy]
shift = 0.0
```

Combined band/DOS label example:

```toml
[band.combined_plot]
dos_width_ratio = 0.30
dos_xlabel = "DOS (states/eV)"  # x label for the rotated DOS panel
dos_ylabel = ""                 # y label for the rotated DOS panel
```

Plot files are opt-in through `[band.plot].targets`; `--no-plot` still disables
all plotting even when targets are present:

```toml
[band.plot]
targets = ["band", "dos", "combined", "fatband", "pdos", "fatband_pdos", "oam", "sam"]
size = [900, 600]
energy_range = [-3.0, 3.0]
font_size = 18
```

Allowed targets are `band`, `dos`, `combined`, `fatband`, `pdos`,
`fatband_pdos`, `oam`, and `sam`. Omit `targets`, or set it to an empty array, to write
data files without plot files. Projection targets require
`[band.projection].enabled = true`; `oam` requires `[band.oam].enabled = true`;
`sam` requires `[band.sam].enabled = true`.

DOS plots include a dotted center-of-mass marker. In standalone DOS plots the
marker is vertical on the energy axis; in combined plots it is horizontal in the
rotated DOS panel. Spin-resolved DOS uses the total `dos_up + dos_down` center,
and PDOS plots use each projection group's color for its own center marker.

Spin-colored band and spin-resolved DOS example:

```toml
[band.spin]
enabled = true
layout = "qe"  # "qe" or "vasp544"
colors = ["#1f77b4", "#d62728"]
```

`layout = "qe"` uses QE/interleaved spin ordering (`up, down, up, down, ...`).
`layout = "vasp544"` uses VASP 5.4.4 block spin ordering
(`all up, then all down`). The layout is explicit; `band` does not infer it
from the Hamiltonian. VASP block-ordered HR, wsvec, and projection indices are
canonicalized to the internal QE/interleaved order before diagonalization.
Band plotting follows the sorted eigenvalue column order for every layout. Spin
coloring and spin-resolved DOS use the canonical internal order
(`up, down, up, down, ...`) after input canonicalization.

Orbital projection with explicit Wannier indices:

```toml
[band.projection]
enabled = true
mode = "index_groups"
plot_style = "colorbar"                 # "colorbar" or "empty_circle"
colorbar_colormap = "viridis"           # any Plots.jl color gradient name
groups = [
  ["Ru_t2g", [1, 2, 3], "#1f77b4"],
  ["Cl_p", [4, 5, 6], "#2ca02c"],
]
```

Orbital projection from `wannier90.win`:

```toml
[band.projection]
enabled = true
mode = "win_groups"
win = "wannier90.win"
groups = [
  ["Mo_sp", ["Mo"], ["s", "p"], "red"],
  ["Mo_t2g", ["Mo"], ["t2g"], "blue"],
]
```

Supported projection modes:

- `index_groups`: groups are direct 1-based Wannier orbital indices.
- `win_groups`: groups are selected from the `begin projections` block in `wannier90.win`.
- `color_group` is optional. If omitted, projected band plots use every group in
  `groups`. Set `color_group = ["label_a"]` only when you want to plot a subset
  of groups.
- `groups` may use compact rows. For `index_groups`, use
  `[label, indices, color]`; for `win_groups`, use
  `[label, species, orbitals, color]`. The `[[band.projection.groups]]`
  table form remains supported for advanced win selectors such as `site_labels`,
  numeric `sites`, or `spin`.

Supported projected-band plot styles:

- When requested by `plot.targets`, projected band plots overlay the selected
  projection groups; PDOS plots always show all groups.
- `plot_style = "colorbar"`: fixed-size markers are colored by each selected group weight. With multiple groups, each group uses its own color gradient.
- `plot_style = "empty_circle"`: hollow circles are drawn with size proportional to each selected group weight and stroke color from the group. Use `circle_max_size` and `circle_stroke_width` to tune the marker geometry.
- `colorbar_colormap = "plasma"` changes the named gradient used by colorbar mode. For a custom gradient, use `colorbar_colors = ["white", "#1f77b4"]`.

Projection data outputs are separate from the existing band/DOS files. Existing
`bands.dat` and `dos.dat` are not changed.
PDOS curves are normalized per atom for each projection group when the projection source contains site metadata (`win_groups`). `index_groups` has no atom metadata, so its atom count is 1.

Projected atomic OAM uses the same projection basis metadata. It is enabled only
for band-path runs and requires `[band.projection]` with `mode = "win_groups"`:

```toml
[band.oam]
enabled = true
orbitals = [
  ["Ru1", "t2g"],
  ["Cl1", "p"],
]
degeneracy_tol = 1.0e-4
plot_components = ["Lz"]
```

Each `orbitals` row is `[site_label, shell]`. The first value is matched only as
a site label; species fallback is not used. The first version supports the shell
tokens `s`, `p`, `d`, `t2g`, and `eg`. `plot_components` defaults to `["Lz"]`
and may include `Lx`, `Ly`, `Lz`, `L_norm`, or `L2`. OAM writes `bands_oam.dat`,
derived from `[band.output] bands_data`. OAM-colored band plots are written only
when `plot.targets` includes `oam`; their paths are derived from
`[band.output] bands_plot`, such as `bands_oam_lz.png`. There are no TOML output
paths for these derived files. The values are projected local atomic OAM in the
chosen Wannier basis, not Berry-phase orbital magnetization. At exact
degeneracies, per-band OAM is gauge-dependent; the first version warns but does
not write degenerate-cluster traces.

Spin angular momentum expectation uses the same projection basis metadata and is
enabled with `[band.sam]`. It requires an explicit `[band.spin] layout` but does
not require `spin.enabled = true`:

```toml
[band.spin]
layout = "vasp544" # or "qe"

[band.sam]
enabled = true
degeneracy_tol = 1.0e-4
plot_components = ["Sz"]
```

SAM writes `bands_sam.dat`, derived from `[band.output] bands_data`, with `Sx`,
`Sy`, `Sz`, `S_norm`, and `S2`. `plot_components` defaults to `["Sz"]` and may
include `Sx`, `Sy`, `Sz`, `S_norm`, or `S2`. SAM-colored band plots are written
only when `plot.targets` includes `sam`; their paths are derived from
`[band.output] bands_plot`, such as `bands_sam_sz.png`. Nonmagnetic or
unpolarized basis metadata is an error because the Pauli blocks require complete
`(site, orbital, up/dn)` pairs.

Default data outputs and plot destinations are resolved relative to the input
TOML directory. Plot destinations are used only when the matching
`plot.targets` entry is present.

| Table/key | Default |
| --- | --- |
| `[band.output] bands_data` | `outputs/data/bands.dat` |
| `[band.output] dos_data` | `outputs/data/dos.dat` |
| `[band.output] bands_plot` | `outputs/plots/bands.png` |
| `[band.output] dos_plot` | `outputs/plots/dos.png` |
| `[band.output] combined_plot` | `outputs/plots/band_dos.png` |
| `[band.projection] weights_data` | `outputs/data/projection/band_projection_weights.dat` |
| `[band.projection] pdos_data` | `outputs/data/projection/pdos.dat` |
| `[band.projection] projected_bands_plot` | `outputs/plots/projection/bands_projected.png` |
| `[band.projection] pdos_plot` | `outputs/plots/projection/pdos.png` |
| `[band.projection] projected_combined_plot` | `outputs/plots/projection/band_pdos.png` |
| `[band.oam] derived data` | `outputs/data/bands_oam.dat` |
| `[band.oam] derived plot` | `outputs/plots/bands_oam_lz.png` |
| `[band.sam] derived data` | `outputs/data/bands_sam.dat` |
| `[band.sam] derived plot` | `outputs/plots/bands_sam_sz.png` |

Run examples:

```bash
julia main.jl band --input input.band.toml
```

```bash
julia main.jl band --input input.band.toml --no-plot
```

SAM example:

```bash
julia spectrum/main.jl spectrum/examples/sam/bi2se3/input.toml --no-plot
```

## Dependencies

Internal dependencies:

- `../shared/CrystalCell.jl`
- `../shared/PoscarIO.jl`
- `../shared/LatticeIO.jl`
- `../shared/WannierHrIO.jl`
- `../shared/WannierWsvecIO.jl`
- `../shared/WannierWsvecGenerate.jl`
- `../shared/WannierKspace.jl`
- `../shared/WannierEigensystem.jl`
- `../shared/OrbitalProjection.jl`
- `../shared/AtomicOam.jl`
- `../shared/AtomicSpin.jl`
- `../shared/Win90Basis.jl`
- `projection/Model.jl`
- `projection/InputGroups.jl`
- `projection/Input.jl`
- `PlotInput.jl`
- `projection/Core.jl`
- `InputIO.jl`
- `KPath.jl`
- `Execution.jl`
- `observables/Oam.jl`
- `observables/Sam.jl`
- `observables/Output.jl`
- `observables/Plot.jl`
- `../shared/SpinLayout.jl`
- `Bands.jl`
- `Dos.jl`
- `Output.jl`
- `PlotService.jl`
- `PlotBands.jl`

External dependencies:

- `Plots`
- Julia `Printf`
- Julia threading support

## Minimal Code Snippet

CLI usage:

```bash
julia main.jl band --input input.band.toml --no-plot
```

Service-layer usage:

```julia
include("spectrum/band/CLI.jl")

cfg = Main.BandCLI.InputIO.read_input("spectrum/examples/graphene/input.toml")
result = Main.BandCLI.Service.run(cfg; make_plot=false)
println(result.hermiticity_ok)
```

## Run And Test

Run the module:

```bash
julia main.jl band --input input.band.toml --no-plot
```

Run tests:

```bash
julia --project=spectrum --startup-file=no spectrum/test/runtests.jl
```

## Common Pitfalls

- `wsvec` corrections are opt-in. Omitting `run.wsvec` and setting `wsvec = ""` both disable them; a non-empty path must exist and cover every `(R, i, j)` entry.
- Cartesian `KPOINTS` require `run.structure` (POSCAR/CONTCAR or `wannier90.win`) so coordinates can be reduced in the input cell's reciprocal basis.
- Reduced `KPOINTS` are interpreted directly in the reciprocal basis of the input cell used by `hr.dat`; HR `R` vectors are not re-based automatically. If `run.structure` is not the same cell convention as `hr.dat`, the K-path distance axis can be distorted.
- Set `run.structure = ""` only when you intentionally want reduced-coordinate
  path distances without a lattice.
- Projection weights are computed from raw eigenvectors. `energy.shift` only changes eigenvalues and DOS energies.
- In `win_groups`, user-requested orbitals are validated against `wannier90.win`; asking for a valid but absent orbital such as `Ru:p` when the win file only has `Ru:d` raises an error instead of producing an empty group.
- At exact degeneracies, individual band projection weights can depend on the eigenvector gauge; group sums over a degenerate subspace are more stable than single-band colors.
- At exact degeneracies, individual band OAM values can also depend on the
  eigenvector gauge. `band.oam` emits per-band values and warnings, not cluster
  traces.
- `band.oam` is a projected local atomic observable in the Wannier basis, not an
  orbital magnetization calculation.
- Unknown `[band.*]` tables and unsupported keys fail early; keep experimental
  notes outside the `[band]` namespace.
- Plot generation depends on `Plots`, so environment issues can break plotting even when computation works.
- `energy.shift = 0.0` means values are raw eigenvalues, not automatically shifted to a Fermi level.

## Reading Order

1. `CLI.jl`
2. `InputIO.jl`
3. `PlotInput.jl`
4. `Service.jl`
5. `KPath.jl`
6. `Bands.jl`
7. `Dos.jl`
8. `projection/Core.jl`
9. `Output.jl`
10. `PlotService.jl`
11. `../shared/SpinLayout.jl`
12. `PlotBands.jl`
13. `observables/Plot.jl`
14. `../shared/WannierEigensystem.jl`
15. `../shared/OrbitalProjection.jl`
16. `../shared/AtomicOam.jl`
17. `../shared/AtomicSpin.jl`
18. `observables/Sam.jl`
19. `../shared/Win90Basis.jl`
