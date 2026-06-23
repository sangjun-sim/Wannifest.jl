module InputIO

using TOML

using ..InputParsing: namespaced_root, optional_string, optional_table, parse_float_pair
using ..InputParsing: parse_int_pair, reject_unknown_keys, reject_unknown_sibling_tables
using ..InputParsing: required_bool, required_table
using ..InputParsing: resolve_path
using ..Model: EnergyConfig, OutputConfig, PlaneConfig, PlotConfig, RunConfig, RunFiles
using ..SpinLayout

export read_input

const DEFAULT_HR = "wannier90_hr.dat"
const DEFAULT_STRUCTURE = "POSCAR"
const DEFAULT_OUTPUT_DIR = "outputs"
const ALLOWED_ROOT_TABLES = Set(("run", "plane", "energy", "spin", "plot"))
const ALLOWED_RUN_KEYS = Set(("hr", "structure", "wsvec", "hermiticity_tol", "verbose"))
const ALLOWED_PLANE_KEYS = Set(("axes", "fixed_axis", "fixed_value", "range_x", "range_y", "mesh"))
const ALLOWED_ENERGY_KEYS = Set(("shift", "bands"))
const ALLOWED_SPIN_KEYS = Set(("layout",))
const ALLOWED_PLOT_KEYS = Set(("mode", "interactive", "size", "energy_range", "colormap", "contour_levels"))
const ALLOWED_PLOT_MODES = Set((:surface, :contour, :heatmap, :both))

function _table_or_empty(root, key::AbstractString)
    tbl = optional_table(root, key)
    return isnothing(tbl) ? Dict{String, Any}() : tbl
end

function _numeric_value(tbl, key::AbstractString, default::Real; context::AbstractString)
    value = get(tbl, key, default)
    value isa Real || error("$context.$key must be numeric")
    return Float64(value)
end

function _integer_value(tbl, key::AbstractString, default::Integer; context::AbstractString)
    value = get(tbl, key, default)
    value isa Integer || error("$context.$key must be an integer")
    return Int(value)
end

function _bool_value(tbl, key::AbstractString, default::Bool; context::AbstractString)
    return haskey(tbl, key) ? required_bool(tbl, key; context=context) : default
end

function _string_value(tbl, key::AbstractString, default::AbstractString; context::AbstractString)
    return optional_string(tbl, key; default=default, context=context)
end

function _resolve_existing_file(base_dir::AbstractString, tbl, key::AbstractString, default::AbstractString)
    raw = _string_value(tbl, key, default; context="contour.run")
    path = resolve_path(base_dir, raw; empty_value="")
    isempty(path) && error("contour.run.$key cannot be empty")
    isfile(path) || error("contour.run.$key points to missing file: $path")
    return path
end

function _resolve_optional_existing_file(base_dir::AbstractString, tbl, key::AbstractString)
    haskey(tbl, key) || return nothing
    raw = optional_string(tbl, key; default="", context="contour.run")
    path = resolve_path(base_dir, raw; empty_value=nothing)
    isnothing(path) && return nothing
    isfile(path) || error("contour.run.$key points to missing file: $path")
    return path
end

function _resolve_structure(base_dir::AbstractString, run_tbl)
    if haskey(run_tbl, "structure")
        raw = optional_string(run_tbl, "structure"; default="", context="contour.run")
        path = resolve_path(base_dir, raw; empty_value=nothing)
        isnothing(path) && return ""
        isfile(path) || error("contour.run.structure points to missing file: $path")
        return path
    end
    default_path = resolve_path(base_dir, DEFAULT_STRUCTURE; empty_value="")
    return isfile(default_path) ? default_path : ""
end

function _axis_index(raw, key::AbstractString)::Int
    raw isa AbstractString || error("$key must be a string")
    label = lowercase(strip(String(raw)))
    label == "kx" && return 1
    label == "ky" && return 2
    label == "kz" && return 3
    error("$key must be one of kx, ky, kz")
end

function _parse_axes(raw)
    raw isa AbstractVector || error("contour.plane.axes must be an array of two axis labels")
    length(raw) == 2 || error("contour.plane.axes must contain exactly two axis labels")
    axes = (_axis_index(raw[1], "contour.plane.axes[1]"), _axis_index(raw[2], "contour.plane.axes[2]"))
    axes[1] != axes[2] || error("contour.plane.axes entries must be distinct")
    return axes
end

function _remaining_axis(x_axis::Int, y_axis::Int)::Int
    for axis in 1:3
        axis != x_axis && axis != y_axis && return axis
    end
    error("contour.plane.axes must leave one fixed axis")
end

function _parse_mesh(raw)::Tuple{Int, Int}
    raw isa AbstractString || error("contour.plane.mesh must be a string like \"101x101\"")
    parts = split(lowercase(strip(String(raw))), 'x')
    length(parts) == 2 || error("contour.plane.mesh must have form NxM")
    mesh = (parse(Int, parts[1]), parse(Int, parts[2]))
    all(>(0), mesh) || error("contour.plane.mesh entries must be positive")
    return mesh
end

function _parse_plane_config(plane_tbl)::PlaneConfig
    reject_unknown_keys(plane_tbl, ALLOWED_PLANE_KEYS, "contour.plane")
    x_axis, y_axis = _parse_axes(get(plane_tbl, "axes", ["kx", "ky"]))
    fixed_axis = haskey(plane_tbl, "fixed_axis") ?
        _axis_index(plane_tbl["fixed_axis"], "contour.plane.fixed_axis") :
        _remaining_axis(x_axis, y_axis)
    fixed_axis != x_axis && fixed_axis != y_axis ||
        error("contour.plane.fixed_axis must be distinct from contour.plane.axes")

    range_x = parse_float_pair(get(plane_tbl, "range_x", [-0.5, 0.5]), "contour.plane.range_x")
    range_y = parse_float_pair(get(plane_tbl, "range_y", [-0.5, 0.5]), "contour.plane.range_y")
    range_x[1] < range_x[2] || error("contour.plane.range_x must have min < max")
    range_y[1] < range_y[2] || error("contour.plane.range_y must have min < max")

    return PlaneConfig(
        x_axis,
        y_axis,
        fixed_axis,
        _numeric_value(plane_tbl, "fixed_value", 0.0; context="contour.plane"),
        range_x,
        range_y,
        _parse_mesh(get(plane_tbl, "mesh", "101x101")),
    )
end

function _parse_bands(raw)::Vector{Int}
    raw isa AbstractVector || error("contour.energy.bands must be an array of positive integers")
    bands = Int[]
    for (i, value) in enumerate(raw)
        value isa Integer || error("contour.energy.bands[$i] must be an integer")
        value > 0 || error("contour.energy.bands[$i] must be positive")
        push!(bands, Int(value))
    end
    isempty(bands) && error("contour.energy.bands cannot be empty")
    return bands
end

function _parse_energy_config(energy_tbl)::EnergyConfig
    reject_unknown_keys(energy_tbl, ALLOWED_ENERGY_KEYS, "contour.energy")
    return EnergyConfig(
        _numeric_value(energy_tbl, "shift", 0.0; context="contour.energy"),
        _parse_bands(get(energy_tbl, "bands", [1])),
    )
end

function _parse_plot_config(plot_tbl)::PlotConfig
    reject_unknown_keys(plot_tbl, ALLOWED_PLOT_KEYS, "contour.plot")
    mode = Symbol(lowercase(_string_value(plot_tbl, "mode", "both"; context="contour.plot")))
    mode in ALLOWED_PLOT_MODES || error("Unsupported contour.plot.mode: $mode")
    size = parse_int_pair(get(plot_tbl, "size", [900, 700]), "contour.plot.size")
    all(>(0), size) || error("contour.plot.size entries must be positive")
    energy_range = parse_float_pair(get(plot_tbl, "energy_range", [-3.0, 3.0]), "contour.plot.energy_range")
    energy_range[1] < energy_range[2] || error("contour.plot.energy_range must have min < max")
    levels = _integer_value(plot_tbl, "contour_levels", 40; context="contour.plot")
    levels > 0 || error("contour.plot.contour_levels must be positive")
    return PlotConfig(
        mode,
        _bool_value(plot_tbl, "interactive", false; context="contour.plot"),
        size,
        energy_range,
        _string_value(plot_tbl, "colormap", "viridis"; context="contour.plot"),
        levels,
    )
end

function _parse_spin_layout(spin_tbl)::Symbol
    reject_unknown_keys(spin_tbl, ALLOWED_SPIN_KEYS, "contour.spin")
    return SpinLayout.parse_layout(get(spin_tbl, "layout", nothing); context="contour.spin.layout")
end

function _resolve_output_dir(base_dir::AbstractString, output_dir_override)
    isnothing(output_dir_override) && return normpath(joinpath(base_dir, DEFAULT_OUTPUT_DIR))
    text = strip(String(output_dir_override))
    isempty(text) && error("--output-dir cannot be empty")
    return resolve_path(base_dir, text; empty_value="")
end

function read_input(path::AbstractString; output_dir_override=nothing)::RunConfig
    cfg = TOML.parsefile(path)
    base_dir = dirname(abspath(path))
    root = namespaced_root(cfg, "contour")
    reject_unknown_sibling_tables(root, ALLOWED_ROOT_TABLES, "contour")

    run_tbl = required_table(root, "run"; context="contour")
    reject_unknown_keys(run_tbl, ALLOWED_RUN_KEYS, "contour.run")
    plane_tbl = required_table(root, "plane"; context="contour")
    energy_tbl = _table_or_empty(root, "energy")
    spin_tbl = _table_or_empty(root, "spin")
    plot_tbl = _table_or_empty(root, "plot")

    files = RunFiles(
        _resolve_existing_file(base_dir, run_tbl, "hr", DEFAULT_HR),
        _resolve_structure(base_dir, run_tbl),
        _resolve_optional_existing_file(base_dir, run_tbl, "wsvec"),
    )

    hermiticity_tol = _numeric_value(run_tbl, "hermiticity_tol", 1.0e-8; context="contour.run")
    hermiticity_tol > 0 || error("contour.run.hermiticity_tol must be positive")
    return RunConfig(
        files,
        _parse_plane_config(plane_tbl),
        _parse_energy_config(energy_tbl),
        _parse_plot_config(plot_tbl),
        OutputConfig(_resolve_output_dir(base_dir, output_dir_override)),
        _parse_spin_layout(spin_tbl),
        hermiticity_tol,
        _bool_value(run_tbl, "verbose", true; context="contour.run"),
    )
end

end
