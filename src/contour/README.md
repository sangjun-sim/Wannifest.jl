# spectrum contour

`spectrum/contour` samples a two-dimensional reduced-coordinate k-plane and
writes `(kx, ky, kz, band, energy)` data plus optional PNG plots.

Run from the repository root:

```bash
julia spectrum/main.jl contour --input path/to/input.contour.toml
```

Skip plot generation:

```bash
julia spectrum/main.jl contour --input path/to/input.contour.toml --no-plot
```

Write outputs under a custom directory:

```bash
julia spectrum/main.jl contour --input path/to/input.contour.toml --output-dir outputs/custom
```

## Input

Example:

```toml
[contour.run]
hr = "wannier90_hr.dat"
structure = ""
wsvec = ""
hermiticity_tol = 1.0e-8
verbose = true

[contour.plane]
axes = ["kx", "ky"]
fixed_axis = "kz"
fixed_value = 0.0
range_x = [-0.5, 0.5]
range_y = [-0.5, 0.5]
mesh = "101x101"

[contour.energy]
shift = 0.0
bands = [1]

[contour.spin]
layout = "qe"

[contour.plot]
mode = "both"
interactive = false
size = [900, 700]
energy_range = [-3.0, 3.0]
colormap = "viridis"
contour_levels = 40
```

Paths are resolved relative to the input TOML. Output file paths are not part of
the TOML schema; use `--output-dir` only when the default directory should move.

`[contour.spin].layout` only describes the source HR/wsvec ordering
(`qe` or `vasp544`). The current contour module visualizes sorted band index
surfaces and does not compute spin expectation maps.

## Outputs

Default outputs are written relative to the input TOML directory:

- `outputs/data/contour_energy_surface.dat`
- `outputs/plots/contour_surface.png`
- `outputs/plots/contour_contour.png`
- `outputs/plots/contour_heatmap.png`

If `interactive = true`, matching Plotly HTML files are also written and opened
in the default browser:

- `outputs/plots/contour_surface.html`
- `outputs/plots/contour_contour.html`
- `outputs/plots/contour_heatmap.html`

The browser plots support zoom, pan, hover, and 3D rotation for surface mode.
`--no-plot` disables both PNG saving and interactive HTML opening. All selected
bands are drawn into the same plot for each plot mode.
