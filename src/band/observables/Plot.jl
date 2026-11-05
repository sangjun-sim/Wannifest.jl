module ObservablesPlot

using Plots

using ..Model: BandResult, RunConfig

export oam_plot_path, plot_oam_bands
export sam_plot_path, plot_sam_bands

const OAM_COMPONENT_INDEX = Dict(
    :Lx => 1,
    :Ly => 2,
    :Lz => 3,
    :L_norm => 4,
    :L2 => 5,
)
const SAM_COMPONENT_INDEX = Dict(
    :Sx => 1,
    :Sy => 2,
    :Sz => 3,
    :S_norm => 4,
    :S2 => 5,
)
const OAM_SIGNED_COMPONENTS = Set((:Lx, :Ly, :Lz))
const SAM_SIGNED_COMPONENTS = Set((:Sx, :Sy, :Sz))

function oam_plot_path(bands_plot_path::AbstractString, component::Symbol=:Lz)
    root, ext = splitext(String(bands_plot_path))
    suffix = lowercase(string(component))
    return string(root, "_oam_", suffix, isempty(ext) ? ".png" : ext)
end

function sam_plot_path(bands_plot_path::AbstractString, component::Symbol=:Sz)
    root, ext = splitext(String(bands_plot_path))
    suffix = lowercase(string(component))
    return string(root, "_sam_", suffix, isempty(ext) ? ".png" : ext)
end

function _subplot_margins(config::RunConfig)
    left, right, top, bottom = config.plot.margins
    return (
        margin=0 * Plots.mm,
        left_margin=left * Plots.mm,
        right_margin=right * Plots.mm,
        top_margin=top * Plots.mm,
        bottom_margin=bottom * Plots.mm,
    )
end

function _font_kwargs(config::RunConfig)
    return (
        guidefontsize=config.plot.font_size,
        tickfontsize=max(config.plot.font_size - 2, 1),
    )
end

# High-symmetry guide lines go in before any data series so the bands stay on
# top of them.
function _plot_tick_vlines!(plt, result::BandResult)
    for xpos in result.tick_positions
        vline!(plt, [xpos]; color=:gray, linewidth=0.8, linestyle=:dash, label=false)
    end
    return plt
end

function _plot_fermi_hline!(plt, config::RunConfig)
    config.plot.fermi_line || return plt
    level = something(config.plot.fermi_level, 0.0)
    lo, hi = config.plot.energy_range
    (isfinite(level) && lo <= level <= hi) || return plt
    hline!(plt, [level]; color=:gray, linewidth=0.9, linestyle=:dash, label=false)
    return plt
end

function _oam_component_values(result::BandResult, component::Symbol)
    isnothing(result.oam) && error("BandResult has no OAM data")
    index = get(OAM_COMPONENT_INDEX, component, nothing)
    isnothing(index) && error("Unsupported OAM plot component: $component")
    return result.oam[:, :, index]
end

function _sam_component_values(result::BandResult, component::Symbol)
    isnothing(result.sam) && error("BandResult has no SAM data")
    index = get(SAM_COMPONENT_INDEX, component, nothing)
    isnothing(index) && error("Unsupported SAM plot component: $component")
    return result.sam[:, :, index]
end

function _component_clims(values, component::Symbol, signed_components)
    if component in signed_components
        limit = maximum(abs, values)
        limit == 0.0 && (limit = 1.0)
        return (-limit, limit)
    end
    limit = maximum(values)
    limit == 0.0 && (limit = 1.0)
    return (0.0, limit)
end

function _plot_observable_bands(
    result::BandResult,
    config::RunConfig,
    values,
    clims;
    path::AbstractString,
    save::Bool,
)
    plt = plot(
        size=config.plot.size,
        legend=false,
        dpi=200,
        framestyle=:box,
        grid=false;
        _font_kwargs(config)...,
        _subplot_margins(config)...,
    )

    _plot_tick_vlines!(plt, result)

    nbands = size(result.eigenvalues, 2)
    colorbar_pending = true
    for ib in 1:nbands
        for seg_range in result.segment_ranges
            plot!(
                plt,
                result.distances[seg_range],
                result.eigenvalues[seg_range, ib];
                color=:gray,
                alpha=0.35,
                linewidth=max(config.band_plot.linewidth * 0.8, 0.5),
            )
            scatter!(
                plt,
                result.distances[seg_range],
                result.eigenvalues[seg_range, ib];
                marker_z=values[seg_range, ib],
                markersize=3,
                markerstrokewidth=0,
                color=cgrad(:viridis),
                clims=clims,
                colorbar=colorbar_pending,
            )
            colorbar_pending = false
        end
    end

    _plot_fermi_hline!(plt, config)
    xlabel!(plt, "")
    ylabel!(plt, config.band_plot.ylabel)
    xticks!(plt, result.tick_positions, result.tick_labels)
    xlims!(plt, extrema(result.distances))
    ylims!(plt, config.plot.energy_range)
    if save
        mkpath(dirname(path))
        savefig(plt, path)
    end
    return plt
end

function plot_oam_bands(result::BandResult, config::RunConfig; component::Symbol=:Lz, save::Bool=true)
    values = _oam_component_values(result, component)
    clims = _component_clims(values, component, OAM_SIGNED_COMPONENTS)
    return _plot_observable_bands(
        result,
        config,
        values,
        clims;
        path=oam_plot_path(config.output.bands_plot, component),
        save=save,
    )
end

function plot_sam_bands(result::BandResult, config::RunConfig; component::Symbol=:Sz, save::Bool=true)
    values = _sam_component_values(result, component)
    clims = _component_clims(values, component, SAM_SIGNED_COMPONENTS)
    return _plot_observable_bands(
        result,
        config,
        values,
        clims;
        path=sam_plot_path(config.output.bands_plot, component),
        save=save,
    )
end

end
