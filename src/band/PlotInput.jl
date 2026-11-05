module PlotInput

using ..InputParsing: integer_value, bool_value
using ..InputParsing: parse_float_pair, parse_int_pair, reject_unknown_keys
using ..Model: PlotConfig

export parse_plot_config, validate_plot_targets

const ALLOWED_PLOT_KEYS =
    Set(("interactive", "size", "energy_range", "font_size", "margins", "targets", "fermi_line", "fermi_level"))
const PLOT_TARGETS = Dict(
    "band" => :band,
    "dos" => :dos,
    "combined" => :combined,
    "fatband" => :fatband,
    "pdos" => :pdos,
    "fatband_pdos" => :fatband_pdos,
    "oam" => :oam,
    "sam" => :sam,
    "compare" => :compare,
)

function _parse_margins(tbl)::NTuple{4, Float64}
    raw = get(tbl, "margins", [6.0, 2.0, 2.0, 5.0])
    raw isa AbstractVector || error("plot.margins must be an array of four values in mm")
    length(raw) == 4 || error("plot.margins must contain [left, right, top, bottom] in mm")
    all(value -> value isa Real, raw) || error("plot.margins must contain only numbers")
    margins = Tuple(Float64.(raw))
    all(>=(0.0), margins) || error("plot.margins values must be non-negative")
    return margins
end

function _parse_fermi_level(plot_tbl)::Union{Nothing, Float64}
    haskey(plot_tbl, "fermi_level") || return nothing
    value = plot_tbl["fermi_level"]
    value isa Real || error("plot.fermi_level must be a number")
    isfinite(value) || error("plot.fermi_level must be finite")
    return Float64(value)
end

function _parse_plot_targets(plot_tbl)::Vector{Symbol}
    haskey(plot_tbl, "targets") || return Symbol[]
    raw = plot_tbl["targets"]
    raw isa AbstractVector || error("plot.targets must be an array of strings")
    targets = Symbol[]
    seen = Set{Symbol}()
    allowed = join(sort!(collect(keys(PLOT_TARGETS))), ", ")
    for (index, item) in enumerate(raw)
        item isa AbstractString ||
            error("plot.targets must contain only strings (bad entry at index $index)")
        key = lowercase(strip(String(item)))
        isempty(key) && error("plot.targets cannot contain empty strings")
        target = get(PLOT_TARGETS, key, nothing)
        isnothing(target) && error("plot.targets entries must be one of $allowed")
        target in seen && error("plot.targets contains duplicate target $target")
        push!(seen, target)
        push!(targets, target)
    end
    return targets
end

function parse_plot_config(plot_tbl)::PlotConfig
    reject_unknown_keys(plot_tbl, ALLOWED_PLOT_KEYS, "plot")
    interactive = bool_value(plot_tbl, "interactive", false; context="plot")
    plot_size = parse_int_pair(get(plot_tbl, "size", [900, 600]), "plot.size")
    all(>(0), plot_size) || error("plot.size entries must be positive")
    energy_range = parse_float_pair(get(plot_tbl, "energy_range", [-3.0, 3.0]), "plot.energy_range")
    energy_range[1] < energy_range[2] || error("plot.energy_range must have min < max")
    font_size = integer_value(plot_tbl, "font_size", 18; context="plot")
    font_size >= 2 || error("plot.font_size must be >= 2")
    fermi_line = bool_value(plot_tbl, "fermi_line", true; context="plot")
    return PlotConfig(
        interactive,
        plot_size,
        energy_range,
        font_size,
        _parse_margins(plot_tbl),
        _parse_plot_targets(plot_tbl),
        fermi_line,
        _parse_fermi_level(plot_tbl),
    )
end

_has_target(targets::Vector{Symbol}, target::Symbol) = target in targets

function validate_plot_targets(
    targets::Vector{Symbol},
    mode::Symbol,
    projection_config,
    oam_config,
    sam_config;
    reference_enabled::Bool=false,
)
    isempty(targets) && return nothing

    if _has_target(targets, :band) && !(mode in (:bands, :all))
        error("plot target \"band\" requires run.mode=\"bands\" or \"all\"")
    end
    if _has_target(targets, :dos) && !(mode in (:dos, :all))
        error("plot target \"dos\" requires run.mode=\"dos\" or \"all\"")
    end
    if _has_target(targets, :combined) && mode != :all
        error("plot target \"combined\" requires run.mode=\"all\"")
    end
    if _has_target(targets, :fatband)
        mode in (:bands, :all) || error("plot target \"fatband\" requires run.mode=\"bands\" or \"all\"")
        projection_config.enabled || error("plot target \"fatband\" requires band.projection.enabled=true")
    end
    if _has_target(targets, :pdos)
        mode in (:dos, :all) || error("plot target \"pdos\" requires run.mode=\"dos\" or \"all\"")
        projection_config.enabled || error("plot target \"pdos\" requires band.projection.enabled=true")
    end
    if _has_target(targets, :fatband_pdos)
        mode == :all || error("plot target \"fatband_pdos\" requires run.mode=\"all\"")
        projection_config.enabled || error("plot target \"fatband_pdos\" requires band.projection.enabled=true")
    end
    if _has_target(targets, :oam)
        mode in (:bands, :all) || error("plot target \"oam\" requires run.mode=\"bands\" or \"all\"")
        oam_config.enabled || error("plot target \"oam\" requires band.oam.enabled=true")
    end
    if _has_target(targets, :sam)
        mode in (:bands, :all) || error("plot target \"sam\" requires run.mode=\"bands\" or \"all\"")
        sam_config.enabled || error("plot target \"sam\" requires band.sam.enabled=true")
    end
    if _has_target(targets, :compare)
        mode in (:bands, :all) || error("plot target \"compare\" requires run.mode=\"bands\" or \"all\"")
        reference_enabled || error("plot target \"compare\" requires band.reference.enabled=true")
    end
    return nothing
end

end
