# flux

Languages: English | [Korean](<README(ko).md>)

`flux` adds directed complex flux terms to an existing Wannier `hr.dat` file and
writes a Plotly 3D HTML view of the seed bonds. It lives entirely inside
`spectrum/` and reuses `spectrum/shared` support modules.

## Usage

```bash
julia spectrum/main.jl flux --input input.flux.toml --output flux_hr.dat --html flux.html --diagnostic flux_diagnostic.tsv
```

Input example:

```toml
[flux.run]
hr = "wannier90_hr.dat"
win = "wannier90.win"
# poscar = "POSCAR"  # use when win is unavailable

[flux.geometry]
search_bounds = [1, 1, 0]
distance_tol = 1.0e-8

[flux.plot]
interactive = true
cell_bounds = [1, 1, 0]
arrow_styles = [
  ["C1:s", 1.0, "#2ca02c"],
  ["C2:s", 1.0, "#d62728"],
]

# Optional for POSCAR fallback when orbital counts differ by species.
[flux.basis]
orbitals_per_atom = [["Mo", 5], ["S", 3]]
# orbitals_per_atom = [5, 3]  # POSCAR species-group order; useful for duplicate labels

[flux.diagnostic]
continuity = true
continuity_tol = 1.0e-10
plaquettes = [
  ["C1-triangle", [
    ["C1", [0, 0, 0]],
    ["C1", [1, 1, 0]],
    ["C1", [1, 0, 0]],
  ]],
]

[[flux.terms]]
term = [
  [2, ["C1", "C1"], [1, 1, 0], [0.0, -0.02]],
  [2, ["C1", "C1"], [1, 0, 0], [0.0, 0.02]],
  [2, ["C1", "C1"], [0, 1, 0], [0.0, 0.02]],
]
```

`nn` is a distance-shell index, not necessarily a physical neighbor label. Integer
endpoints are 1-based Wannier indices. String endpoints first match `site_label`
values such as `Mo1`; if no site matches, they match species such as `Mo`.
Each term row is `[nn, [from, to], [rx, ry, rz], [re, im]]`, so the seed cell
translation `R` is explicit. Hermitian-equivalent seeds selected by a broad
selector are deduplicated; duplicate explicit rows for the same physical pair
are rejected.

Output paths are CLI options, not TOML keys. The written HR is normalized with
`ndegen = 1`; the HTML shows seed flux bonds and their incoming periodic images,
not Hermitian partners. A row maps to the Wannier matrix element
`<from,0|H|to,R>`, so the positive phase direction is `to@R -> from@0`.
Direction cones are drawn at bond midpoints with enlarged arrowheads so bond
orientation remains visible in the 3D view. Negative imaginary phase values
reverse the visual arrow to `from@0 -> to@R`.
`[flux.plot] cell_bounds`
controls the displayed supercell image range for both atom markers and repeated
periodic seed arrows; when omitted, it falls back to `[flux.geometry] search_bounds`.
`arrow_styles` rows are
`[from_atom_orbital_label, arrow_size, color]` overrides for directed seed
arrows and their incoming periodic images.

`[flux.diagnostic]` writes a TSV report under `outputs/diagnostic/` by default.
Use `--diagnostic path.tsv` to override it or `--no-diagnostic` to suppress it.
Plaquette diagnostics sum only `imag(value)` from matched seed flux edges; they
do not compute the final Hamiltonian Wilson-loop phase. `continuity = true`
also reports per-orbital flow in/out residuals using the same imaginary-sign
direction convention as the arrows.

If `win` is absent, set `poscar`. POSCAR fallback has no projection metadata, so
it assigns HR indices to atoms in POSCAR order. Without `[flux.basis]`, this is
only allowed when `num_wann` divides evenly by the number of atoms. Duplicate
species rows with the same count are accepted. If duplicate POSCAR species
groups need different counts, use the integer-vector form in POSCAR group order.
VASP4-style POSCARs without element names use synthetic species `Type1`,
`Type2`, with site labels like `Type1_1`.
