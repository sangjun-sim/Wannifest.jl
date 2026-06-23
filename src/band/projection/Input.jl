module ProjectionInput

using ..InputParsing: optional_string, required_bool, required_string_vector
using ..InputParsing: reject_unknown_keys, resolve_path
using ..LocalAxisRotation
using ..ProjectionModel:
    ProjectionBasisRotationConfig,
    ProjectionConfig,
    disabled_projection_config
using ..ProjectionInputGroups

export parse_projection_config

const DEFAULT_PROJECTION_WEIGHTS_DATA = "outputs/data/projection/band_projection_weights.dat"
const DEFAULT_PDOS_DATA = "outputs/data/projection/pdos.dat"
const DEFAULT_PROJECTED_BANDS_PLOT = "outputs/plots/projection/bands_projected.png"
const DEFAULT_PDOS_PLOT = "outputs/plots/projection/pdos.png"
const DEFAULT_PROJECTED_COMBINED_PLOT = "outputs/plots/projection/band_pdos.png"
const ALLOWED_PROJECTION_MODES = Set((:index_groups, :win_groups))
const ALLOWED_PROJECTION_PLOT_STYLES = Set((:colorbar, :empty_circle))
const ALLOWED_PROJECTION_KEYS = Set((
    "enabled",
    "mode",
    "color_group",
    "plot_style",
    "colorbar_colormap",
    "colorbar_colors",
    "circle_max_size",
    "circle_stroke_width",
    "weights_data",
    "projected_bands_plot",
    "pdos_data",
    "pdos_plot",
    "projected_combined_plot",
    "win",
    "groups",
    "basis_rotation",
))
const ALLOWED_BASIS_ROTATION_KEYS = Set(("enabled", "local_axes", "strict_t2g", "leakage_tol"))

function _numeric_value(tbl, key::AbstractString, default::Real; context::AbstractString)
    value = get(tbl, key, default)
    value isa Real || error("$context.$key must be numeric")
    return Float64(value)
end

function _bool_value(tbl, key::AbstractString, default::Bool; context::AbstractString)
    return haskey(tbl, key) ? required_bool(tbl, key; context=context) : default
end

function _string_value(tbl, key::AbstractString, default::AbstractString; context::AbstractString)
    return optional_string(tbl, key; default=default, context=context)
end

function _optional_string_vector(tbl, key::AbstractString; context::AbstractString)
    haskey(tbl, key) || return String[]
    return required_string_vector(tbl, key; context=context)
end

function _positive_float(tbl, key::AbstractString, default::Real; context::AbstractString)
    value = _numeric_value(tbl, key, default; context=context)
    value > 0 || error("$context.$key must be positive")
    return value
end

function _parse_projection_plot_style(raw)::Symbol
    raw isa AbstractString || error("projection.plot_style must be a string")
    style = Symbol(strip(String(raw)))
    style in ALLOWED_PROJECTION_PLOT_STYLES ||
        error("Unsupported projection.plot_style: $raw (allowed: colorbar, empty_circle)")
    return style
end

function _projection_output_path(base_dir::AbstractString, tbl, key::AbstractString, default::AbstractString)
    raw = optional_string(tbl, key; default=default, context="projection")
    return resolve_path(base_dir, raw; empty_value="")
end

function _projection_existing_path(
    base_dir::AbstractString,
    tbl,
    key::AbstractString;
    required::Bool,
)
    if !haskey(tbl, key)
        required && error("projection.$key is required")
        return nothing
    end
    raw = optional_string(tbl, key; default="", context="projection")
    path = resolve_path(base_dir, raw; empty_value=nothing)
    if isnothing(path)
        required && error("projection.$key cannot be empty")
        return nothing
    end
    isfile(path) || error("projection.$key points to missing file: $path")
    return path
end

function _parse_projection_basis_rotation(
    proj_tbl,
    mode::Symbol,
)::ProjectionBasisRotationConfig
    haskey(proj_tbl, "basis_rotation") ||
        return ProjectionBasisRotationConfig(false, LocalAxisRotation.AxisSpec[], false, 1.0e-8)
    tbl = proj_tbl["basis_rotation"]
    tbl isa AbstractDict || error("projection.basis_rotation must be a table")
    reject_unknown_keys(tbl, ALLOWED_BASIS_ROTATION_KEYS, "projection.basis_rotation")

    enabled = _bool_value(tbl, "enabled", true; context="projection.basis_rotation")
    strict_t2g = _bool_value(tbl, "strict_t2g", false; context="projection.basis_rotation")
    leakage_tol = _positive_float(tbl, "leakage_tol", 1.0e-8; context="projection.basis_rotation")
    if !enabled
        return ProjectionBasisRotationConfig(false, LocalAxisRotation.AxisSpec[], strict_t2g, leakage_tol)
    end
    mode == :index_groups &&
        error("projection.basis_rotation requires mode=\"win_groups\"")
    haskey(tbl, "local_axes") || error("projection.basis_rotation.local_axes is required when enabled=true")
    local_axes = LocalAxisRotation.parse_axes(tbl["local_axes"])
    return ProjectionBasisRotationConfig(true, local_axes, strict_t2g, leakage_tol)
end

function parse_projection_config(proj_tbl, base_dir::AbstractString)::ProjectionConfig
    reject_unknown_keys(proj_tbl, ALLOWED_PROJECTION_KEYS, "projection")
    enabled = _bool_value(proj_tbl, "enabled", false; context="projection")
    enabled || return disabled_projection_config()

    mode_raw = _string_value(proj_tbl, "mode", "index_groups"; context="projection")
    mode = Symbol(mode_raw)
    mode in ALLOWED_PROJECTION_MODES ||
        error("Unsupported projection.mode: $mode_raw (allowed: index_groups, win_groups)")

    win_path = _projection_existing_path(base_dir, proj_tbl, "win"; required=mode == :win_groups)
    mode == :index_groups && !isnothing(win_path) &&
        error("projection.win is only valid for mode=\"win_groups\"")
    basis_rotation = _parse_projection_basis_rotation(proj_tbl, mode)

    groups = ProjectionInputGroups.parse_projection_groups(proj_tbl, mode)
    labels = [group.label for group in groups]

    color_group = ProjectionInputGroups.projection_color_groups(proj_tbl, labels)
    plot_style = _parse_projection_plot_style(get(proj_tbl, "plot_style", "colorbar"))
    colorbar_colormap = _string_value(proj_tbl, "colorbar_colormap", "viridis"; context="projection")
    colorbar_colors = _optional_string_vector(proj_tbl, "colorbar_colors"; context="projection")
    circle_max_size = _positive_float(proj_tbl, "circle_max_size", 9.0; context="projection")
    circle_stroke_width = _positive_float(proj_tbl, "circle_stroke_width", 1.0; context="projection")

    return ProjectionConfig(
        true,
        mode,
        color_group,
        plot_style,
        colorbar_colormap,
        colorbar_colors,
        circle_max_size,
        circle_stroke_width,
        _projection_output_path(base_dir, proj_tbl, "weights_data", DEFAULT_PROJECTION_WEIGHTS_DATA),
        _projection_output_path(base_dir, proj_tbl, "projected_bands_plot", DEFAULT_PROJECTED_BANDS_PLOT),
        _projection_output_path(base_dir, proj_tbl, "pdos_data", DEFAULT_PDOS_DATA),
        _projection_output_path(base_dir, proj_tbl, "pdos_plot", DEFAULT_PDOS_PLOT),
        _projection_output_path(base_dir, proj_tbl, "projected_combined_plot", DEFAULT_PROJECTED_COMBINED_PLOT),
        win_path,
        groups,
        basis_rotation,
    )
end

end
