module InputIO

using TOML

using ..InputParsing: namespaced_root, reject_unknown_keys, reject_unknown_sibling_tables
using ..InputParsing: required_table, optional_table, required_string, optional_string
using ..InputParsing: required_bool, resolve_path
using ..Model: ArrowStyle, FluxConfig, FluxEndpoint, FluxRow, FluxTerm, GeometryConfig, PlotConfig
using ..Model: BasisConfig, DiagnosticConfig, DiagnosticVertex, PlaquetteDiagnostic, RunFiles, SpinConfig
using ..SpinLayout

export read_input

const ALLOWED_ROOT_TABLES = Set(("run", "geometry", "plot", "spin", "basis", "diagnostic", "terms"))
const ALLOWED_RUN_KEYS = Set(("hr", "win", "poscar"))
const ALLOWED_GEOMETRY_KEYS = Set(("search_bounds", "distance_tol"))
const ALLOWED_PLOT_KEYS = Set(("interactive", "cell_bounds", "arrow_styles"))
const ALLOWED_SPIN_KEYS = Set(("layout",))
const ALLOWED_BASIS_KEYS = Set(("orbitals_per_atom",))
const ALLOWED_DIAGNOSTIC_KEYS = Set(("enabled", "continuity", "continuity_tol", "plaquettes"))
const ALLOWED_TERM_KEYS = Set(("term",))
const DEFAULT_DISTANCE_TOL = 1.0e-8
const DEFAULT_SEARCH_BOUNDS = (1, 1, 0)

function _resolve_existing_file(base_dir::AbstractString, tbl, key::AbstractString; context::AbstractString)
    raw = required_string(tbl, key; context=context)
    path = resolve_path(base_dir, raw; empty_value="")
    isfile(path) || error("$context.$key points to missing file: $path")
    return path
end

function _resolve_optional_existing_file(base_dir::AbstractString, tbl, key::AbstractString; context::AbstractString)
    raw = optional_string(tbl, key; default=nothing, context=context)
    path = resolve_path(base_dir, raw; empty_value=nothing)
    isnothing(path) && return nothing
    isfile(path) || error("$context.$key points to missing file: $path")
    return String(path)
end

function _parse_vec3_int(raw, key::AbstractString)
    raw isa AbstractVector || error("$key must be a 3-element integer array")
    length(raw) == 3 || error("$key must have length 3")
    all(x -> x isa Integer && !(x isa Bool), raw) || error("$key entries must be integers")
    return (Int(raw[1]), Int(raw[2]), Int(raw[3]))
end

function _parse_endpoint(raw, context::AbstractString)::FluxEndpoint
    raw isa Integer && !(raw isa Bool) && return FluxEndpoint(Int(raw))
    raw isa AbstractString && !isempty(strip(String(raw))) && return FluxEndpoint(String(strip(String(raw))))
    error("$context endpoint must be an integer Wannier index or non-empty string site label")
end

function _parse_complex_pair(raw, context::AbstractString)::ComplexF64
    raw isa AbstractVector || error("$context value must be [re, im]")
    length(raw) == 2 || error("$context value must have length 2")
    raw[1] isa Real && raw[2] isa Real || error("$context value entries must be numeric")
    return ComplexF64(Float64(raw[1]), Float64(raw[2]))
end

function _parse_flux_row(raw, index::Integer)::FluxRow
    context = "flux.terms row $index"
    raw isa AbstractVector || error("$context must be [nn, [from, to], [rx, ry, rz], [re, im]]")
    length(raw) == 4 || error("$context must have length 4")
    raw[1] isa Integer && !(raw[1] isa Bool) || error("$context nn must be an integer")
    nn = Int(raw[1])
    nn >= 1 || error("$context nn must be >= 1")

    endpoints = raw[2]
    endpoints isa AbstractVector || error("$context endpoints must be [from, to]")
    length(endpoints) == 2 || error("$context endpoints must have length 2")
    from = _parse_endpoint(endpoints[1], "$context from")
    to = _parse_endpoint(endpoints[2], "$context to")
    typeof(from.value) == typeof(to.value) ||
        error("$context cannot mix integer and string endpoints")
    R = _parse_vec3_int(raw[3], "$context R")
    value = _parse_complex_pair(raw[4], context)
    return FluxRow(nn, from, to, R, value)
end

function _parse_terms(raw)::Vector{FluxTerm}
    raw isa AbstractVector || error("flux.terms must be an array of tables")
    isempty(raw) && error("flux.terms cannot be empty")
    terms = FluxTerm[]
    for (term_index, tbl) in enumerate(raw)
        tbl isa AbstractDict || error("flux.terms[$term_index] must be a table")
        reject_unknown_keys(tbl, ALLOWED_TERM_KEYS, "flux.terms[$term_index]")
        rows_raw = get(tbl, "term", nothing)
        rows_raw isa AbstractVector || error("flux.terms[$term_index].term must be an array")
        isempty(rows_raw) && error("flux.terms[$term_index].term cannot be empty")
        rows = [_parse_flux_row(row, row_index) for (row_index, row) in enumerate(rows_raw)]
        push!(terms, FluxTerm(rows))
    end
    return terms
end

function _parse_geometry(tbl)::GeometryConfig
    reject_unknown_keys(tbl, ALLOWED_GEOMETRY_KEYS, "flux.geometry")
    search_bounds = haskey(tbl, "search_bounds") ?
        _parse_vec3_int(tbl["search_bounds"], "flux.geometry.search_bounds") :
        DEFAULT_SEARCH_BOUNDS
    all(>=(0), search_bounds) || error("flux.geometry.search_bounds entries must be >= 0")
    raw_tol = get(tbl, "distance_tol", DEFAULT_DISTANCE_TOL)
    raw_tol isa Real && !(raw_tol isa Bool) || error("flux.geometry.distance_tol must be numeric")
    distance_tol = Float64(raw_tol)
    distance_tol > 0 || error("flux.geometry.distance_tol must be positive")
    return GeometryConfig(search_bounds, distance_tol)
end

function _parse_positive_float(raw, context::AbstractString)::Float64
    raw isa Real && !(raw isa Bool) || error("$context must be numeric")
    value = Float64(raw)
    value > 0 || error("$context must be positive")
    return value
end

function _parse_nonempty_string(raw, context::AbstractString)::String
    raw isa AbstractString || error("$context must be a string")
    value = strip(String(raw))
    isempty(value) && error("$context cannot be empty")
    return value
end

function _parse_arrow_style(row, index::Integer)::ArrowStyle
    context = "flux.plot.arrow_styles row $index"
    row isa AbstractVector || error("$context must be [atom_orbital, size, color]")
    length(row) == 3 || error("$context must have length 3")
    selector = _parse_nonempty_string(row[1], "$context atom_orbital")
    size = _parse_positive_float(row[2], "$context size")
    color = _parse_nonempty_string(row[3], "$context color")
    return ArrowStyle(selector, size, color)
end

function _parse_arrow_styles(raw)::Vector{ArrowStyle}
    isnothing(raw) && return ArrowStyle[]
    raw isa AbstractVector || error("flux.plot.arrow_styles must be an array of rows")
    return [_parse_arrow_style(row, index) for (index, row) in enumerate(raw)]
end

function _parse_plot(tbl)::PlotConfig
    reject_unknown_keys(tbl, ALLOWED_PLOT_KEYS, "flux.plot")
    interactive = haskey(tbl, "interactive") ? required_bool(tbl, "interactive"; context="flux.plot") : true
    cell_bounds = haskey(tbl, "cell_bounds") ?
        _parse_vec3_int(tbl["cell_bounds"], "flux.plot.cell_bounds") :
        nothing
    isnothing(cell_bounds) || all(>=(0), cell_bounds) ||
        error("flux.plot.cell_bounds entries must be >= 0")
    arrow_styles = _parse_arrow_styles(get(tbl, "arrow_styles", nothing))
    return PlotConfig(interactive, cell_bounds, arrow_styles)
end

function _parse_spin(tbl)::SpinConfig
    reject_unknown_keys(tbl, ALLOWED_SPIN_KEYS, "flux.spin")
    layout = SpinLayout.parse_layout(optional_string(tbl, "layout"; default=nothing, context="flux.spin"))
    return SpinConfig(layout)
end

function _parse_count(raw, context::AbstractString)::Int
    raw isa Integer && !(raw isa Bool) || error("$context count must be an integer")
    count = Int(raw)
    count > 0 || error("$context count must be positive")
    return count
end

function _parse_orbitals_per_atom(raw)::Tuple{Dict{String, Int}, Vector{Int}}
    counts = Dict{String, Int}()
    isnothing(raw) && return counts, Int[]
    if raw isa Integer && !(raw isa Bool)
        Int(raw) > 0 || error("flux.basis.orbitals_per_atom must be positive")
        counts["*"] = Int(raw)
        return counts, Int[]
    end
    raw isa AbstractVector || error("flux.basis.orbitals_per_atom must be an integer or row array")
    if all(item -> item isa Integer && !(item isa Bool), raw)
        group_counts = [_parse_count(item, "flux.basis.orbitals_per_atom group $index")
                        for (index, item) in enumerate(raw)]
        return counts, group_counts
    end
    for (index, row) in enumerate(raw)
        row isa AbstractVector && length(row) == 2 ||
            error("flux.basis.orbitals_per_atom row $index must be [species, count]")
        row[1] isa AbstractString || error("flux.basis.orbitals_per_atom row $index species must be a string")
        species = strip(String(row[1]))
        isempty(species) && error("flux.basis.orbitals_per_atom row $index species cannot be empty")
        count = _parse_count(row[2], "flux.basis.orbitals_per_atom row $index")
        if haskey(counts, species)
            counts[species] == count || error(
                "Conflicting flux.basis.orbitals_per_atom entries for $species; " *
                "use an integer vector like [5, 3] for duplicate POSCAR species groups",
            )
            continue
        end
        counts[species] = count
    end
    return counts, Int[]
end

function _parse_basis(tbl)::BasisConfig
    reject_unknown_keys(tbl, ALLOWED_BASIS_KEYS, "flux.basis")
    counts, group_counts = _parse_orbitals_per_atom(get(tbl, "orbitals_per_atom", nothing))
    return BasisConfig(counts, group_counts)
end

function _parse_diagnostic_vertex(raw, context::AbstractString)::DiagnosticVertex
    raw isa AbstractVector || error("$context must be [orbital_selector, [rx, ry, rz]]")
    length(raw) == 2 || error("$context must have length 2")
    return DiagnosticVertex(_parse_endpoint(raw[1], context), _parse_vec3_int(raw[2], "$context cell"))
end

function _parse_plaquette(row, index::Integer)::PlaquetteDiagnostic
    context = "flux.diagnostic.plaquettes[$index]"
    row isa AbstractVector || error("$context must be [name, vertices]")
    length(row) == 2 || error("$context must have length 2")
    name = _parse_nonempty_string(row[1], "$context name")
    vertices_raw = row[2]
    vertices_raw isa AbstractVector || error("$context vertices must be an array")
    length(vertices_raw) >= 3 || error("$context requires at least 3 vertices")
    vertices = [
        _parse_diagnostic_vertex(vertex, "$context vertex $vertex_index")
        for (vertex_index, vertex) in enumerate(vertices_raw)
    ]
    return PlaquetteDiagnostic(name, vertices)
end

function _parse_plaquettes(raw)::Vector{PlaquetteDiagnostic}
    isnothing(raw) && return PlaquetteDiagnostic[]
    raw isa AbstractVector || error("flux.diagnostic.plaquettes must be an array")
    return [_parse_plaquette(row, index) for (index, row) in enumerate(raw)]
end

function _parse_diagnostic(tbl)::DiagnosticConfig
    reject_unknown_keys(tbl, ALLOWED_DIAGNOSTIC_KEYS, "flux.diagnostic")
    plaquettes = _parse_plaquettes(get(tbl, "plaquettes", nothing))
    continuity = haskey(tbl, "continuity") ?
        required_bool(tbl, "continuity"; context="flux.diagnostic") :
        false
    raw_tol = get(tbl, "continuity_tol", 1.0e-10)
    raw_tol isa Real && !(raw_tol isa Bool) || error("flux.diagnostic.continuity_tol must be numeric")
    continuity_tol = Float64(raw_tol)
    continuity_tol > 0 || error("flux.diagnostic.continuity_tol must be positive")
    default_enabled = continuity || !isempty(plaquettes)
    enabled = haskey(tbl, "enabled") ?
        required_bool(tbl, "enabled"; context="flux.diagnostic") :
        default_enabled
    enabled && !continuity && isempty(plaquettes) &&
        error("flux.diagnostic enabled but no plaquettes or continuity check were requested")
    return DiagnosticConfig(enabled, continuity, continuity_tol, plaquettes)
end

function read_input(path::AbstractString)::FluxConfig
    cfg = TOML.parsefile(path)
    base_dir = dirname(abspath(path))
    root = namespaced_root(cfg, "flux")
    reject_unknown_sibling_tables(root, ALLOWED_ROOT_TABLES, "flux")

    run_tbl = required_table(root, "run"; context="flux")
    reject_unknown_keys(run_tbl, ALLOWED_RUN_KEYS, "flux.run")
    hr_path = _resolve_existing_file(base_dir, run_tbl, "hr"; context="flux.run")
    win_path = _resolve_optional_existing_file(base_dir, run_tbl, "win"; context="flux.run")
    poscar_path = _resolve_optional_existing_file(base_dir, run_tbl, "poscar"; context="flux.run")
    isnothing(win_path) && isnothing(poscar_path) &&
        error("flux.run requires either \"win\" or \"poscar\"")
    files = RunFiles(
        hr_path,
        win_path,
        poscar_path,
    )

    geometry_tbl = optional_table(root, "geometry")
    plot_tbl = optional_table(root, "plot")
    spin_tbl = optional_table(root, "spin")
    basis_tbl = optional_table(root, "basis")
    diagnostic_tbl = optional_table(root, "diagnostic")

    geometry = _parse_geometry(isnothing(geometry_tbl) ? Dict{String, Any}() : geometry_tbl)
    plot = _parse_plot(isnothing(plot_tbl) ? Dict{String, Any}() : plot_tbl)
    spin = _parse_spin(isnothing(spin_tbl) ? Dict{String, Any}() : spin_tbl)
    basis = _parse_basis(isnothing(basis_tbl) ? Dict{String, Any}() : basis_tbl)
    diagnostic = _parse_diagnostic(isnothing(diagnostic_tbl) ? Dict{String, Any}() : diagnostic_tbl)
    terms = _parse_terms(get(root, "terms", nothing))

    return FluxConfig(files, geometry, plot, spin, basis, diagnostic, terms)
end

end
