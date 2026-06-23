module PlotContour

using Plots

using ..Model: EnergySurfaceResult, RunConfig
using ..Output
using ..PlotInteractive

export plot_surface, plot_contour, plot_heatmap, open_interactive

_axis_label(axis::Int) = axis == 1 ? "kx" : axis == 2 ? "ky" : "kz"

function _xy_labels(config::RunConfig)
    return _axis_label(config.plane.x_axis), _axis_label(config.plane.y_axis)
end

function _color(config::RunConfig)
    return Symbol(config.plot.colormap)
end

function _save_plot(plt, config::RunConfig, kind::Symbol, save::Bool)
    save || return plt
    path = Output.default_plot_path(config, kind)
    mkpath(dirname(path))
    savefig(plt, path)
    return plt
end

_band_label(result::EnergySurfaceResult, slot::Int) = "band $(result.bands[slot])"

function _grid_layout(n::Int)
    cols = ceil(Int, sqrt(n))
    rows = ceil(Int, n / cols)
    return (rows, cols)
end

function plot_surface(result::EnergySurfaceResult, config::RunConfig; save::Bool=true)
    xlabel, ylabel = _xy_labels(config)
    z = result.energies[:, :, 1]
    plt = surface(
        result.x_axis,
        result.y_axis,
        z;
        xlabel=xlabel,
        ylabel=ylabel,
        zlabel="E (eV)",
        color=_color(config),
        size=config.plot.size,
        zlims=config.plot.energy_range,
        camera=(35, 30),
        legend=false,
        dpi=200,
    )
    for band_slot in 2:length(result.bands)
        surface!(
            plt,
            result.x_axis,
            result.y_axis,
            result.energies[:, :, band_slot];
            color=_color(config),
            alpha=0.65,
        )
    end
    return _save_plot(plt, config, :surface, save)
end

function plot_contour(result::EnergySurfaceResult, config::RunConfig; save::Bool=true)
    xlabel, ylabel = _xy_labels(config)
    plt = plot(; xlabel=xlabel, ylabel=ylabel, size=config.plot.size, legend=:topright, dpi=200)
    for band_slot in eachindex(result.bands)
        contour!(
            plt,
            result.x_axis,
            result.y_axis,
            result.energies[:, :, band_slot];
            levels=config.plot.contour_levels,
            label=_band_label(result, band_slot),
        )
    end
    return _save_plot(plt, config, :contour, save)
end

function plot_heatmap(result::EnergySurfaceResult, config::RunConfig; save::Bool=true)
    xlabel, ylabel = _xy_labels(config)
    panels = Any[]
    for band_slot in eachindex(result.bands)
        push!(
            panels,
            heatmap(
                result.x_axis,
                result.y_axis,
                result.energies[:, :, band_slot];
                xlabel=xlabel,
                ylabel=ylabel,
                title=_band_label(result, band_slot),
                color=_color(config),
                clims=config.plot.energy_range,
                colorbar=band_slot == length(result.bands),
                dpi=200,
            ),
        )
    end
    plt = plot(panels...; layout=_grid_layout(length(panels)), size=config.plot.size)
    return _save_plot(plt, config, :heatmap, save)
end

function open_interactive(args...; kwargs...)
    return PlotInteractive.open_interactive(args...; kwargs...)
end

end
