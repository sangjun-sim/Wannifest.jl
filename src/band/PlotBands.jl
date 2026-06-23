module PlotBands

using Plots

using ..Model: BandResult, DosResult, RunConfig
using ..SpinLayout

export plot_bands, plot_dos, plot_combined
export plot_projected_bands, plot_projected_dos, plot_projected_combined

function _subplot_margins()
    return (
        margin=2 * Plots.mm,
        left_margin=6 * Plots.mm,
        right_margin=2 * Plots.mm,
        top_margin=2 * Plots.mm,
        bottom_margin=5 * Plots.mm,
    )
end

function _band_color(config::RunConfig, ib::Integer, nbands::Integer)
    return SpinLayout.band_color(config, ib, nbands)
end

function _band_plot_order(config::RunConfig, nbands::Integer)
    return SpinLayout.band_plot_order(config, nbands)
end

function _combined_widths(config::RunConfig)
    dos_width = config.combined_plot.dos_width_ratio / (1 + config.combined_plot.dos_width_ratio)
    return [1 - dos_width, dos_width]
end

_hidden_ticks() = (Float64[], String[])

function _font_kwargs(config::RunConfig)
    return (
        guidefontsize=config.plot.font_size,
        tickfontsize=max(config.plot.font_size - 2, 1),
    )
end

_plot_grid_kwargs() = (grid=false,)

function _center_visible(center::Real, config::RunConfig)
    lo, hi = config.plot.energy_range
    value = Float64(center)
    return isfinite(value) && lo <= value <= hi
end

function _center_marker_kwargs(config::RunConfig; color=:gray)
    return (
        color=color,
        linewidth=max(config.dos_plot.linewidth * 0.55, 0.8),
        linestyle=:dot,
        label=false,
    )
end

function _plot_center_vline!(plt, center::Real, config::RunConfig; color=:gray)
    _center_visible(center, config) || return plt
    vline!(plt, [Float64(center)]; _center_marker_kwargs(config; color=color)...)
    return plt
end

function _plot_center_hline!(plt, center::Real, config::RunConfig; color=:gray)
    _center_visible(center, config) || return plt
    hline!(plt, [Float64(center)]; _center_marker_kwargs(config; color=color)...)
    return plt
end

function plot_bands(result::BandResult, config::RunConfig; save::Bool=true)
    plt = plot(
        size=config.plot.size,
        legend=false,
        dpi=200,
        framestyle=:box;
        _plot_grid_kwargs()...,
        _font_kwargs(config)...,
        _subplot_margins()...,
    )

    nbands = size(result.eigenvalues, 2)
    for ib in _band_plot_order(config, nbands)
        for seg_range in result.segment_ranges
            plot!(
                plt,
                result.distances[seg_range],
                result.eigenvalues[seg_range, ib];
                color=_band_color(config, ib, nbands),
                linewidth=config.band_plot.linewidth,
            )
        end
    end

    for xpos in result.tick_positions
        vline!(plt, [xpos]; color=:gray, linewidth=0.8, linestyle=:dash)
    end
    xlabel!(plt, "")
    ylabel!(plt, config.band_plot.ylabel)
    xticks!(plt, result.tick_positions, result.tick_labels)
    xlims!(plt, extrema(result.distances))
    ylims!(plt, config.plot.energy_range)
    if save
        mkpath(dirname(config.output.bands_plot))
        savefig(plt, config.output.bands_plot)
    end
    return plt
end

function _projection_group_indices(labels::Vector{String}, color_groups::Vector{String})
    isempty(color_groups) && error("color_group cannot be empty")
    length(unique(color_groups)) == length(color_groups) || error("color_group contains duplicate labels")
    indices = Int[]
    missing = String[]
    for color_group in color_groups
        index = findfirst(==(color_group), labels)
        if isnothing(index)
            push!(missing, color_group)
        else
            push!(indices, index)
        end
    end
    isempty(missing) || error("color_group label(s) not in projection labels: $(join(missing, ", "))")
    return indices
end

function _projection_color_gradient(config::RunConfig)
    colors = config.projection.colorbar_colors
    return isempty(colors) ? Symbol(config.projection.colorbar_colormap) : cgrad(colors)
end

function _projection_color_gradient(config::RunConfig, group_color::AbstractString, ngroups::Integer)
    ngroups == 1 && return _projection_color_gradient(config)
    colors = config.projection.colorbar_colors
    return isempty(colors) ? cgrad(["white", String(group_color)]) : cgrad(colors)
end

function _projection_marker_sizes(weights, config::RunConfig)
    return config.projection.circle_max_size .* clamp.(Float64.(weights), 0.0, 1.0)
end

function plot_projected_bands(result::BandResult, config::RunConfig; save::Bool=true)
    projection = result.projection
    isnothing(projection) && error("BandResult has no projection data")
    group_indices = _projection_group_indices(projection.labels, config.projection.color_group)

    plt = plot(
        size=config.plot.size,
        legend=false,
        dpi=200,
        framestyle=:box;
        _plot_grid_kwargs()...,
        _font_kwargs(config)...,
        _subplot_margins()...,
    )

    nbands = size(result.eigenvalues, 2)
    for ib in _band_plot_order(config, nbands)
        for seg_range in result.segment_ranges
            plot!(
                plt,
                result.distances[seg_range],
                result.eigenvalues[seg_range, ib];
                color=:gray,
                alpha=0.35,
                linewidth=max(config.band_plot.linewidth * 0.8, 0.5),
            )
            for (igroup, ig) in enumerate(group_indices)
                weights = projection.weights[seg_range, ib, ig]
                if config.projection.plot_style == :colorbar
                    scatter!(
                        plt,
                        result.distances[seg_range],
                        result.eigenvalues[seg_range, ib];
                        marker_z=weights,
                        markersize=3,
                        markerstrokewidth=0,
                        color=_projection_color_gradient(config, projection.colors[ig], length(group_indices)),
                        clims=(0.0, 1.0),
                        colorbar=igroup == 1,
                    )
                elseif config.projection.plot_style == :empty_circle
                    scatter!(
                        plt,
                        result.distances[seg_range],
                        result.eigenvalues[seg_range, ib];
                        markershape=:circle,
                        markersize=_projection_marker_sizes(weights, config),
                        markercolor=:white,
                        markeralpha=0.3,
                        markerstrokecolor=projection.colors[ig],
                        markerstrokealpha=1.0,
                        markerstrokewidth=config.projection.circle_stroke_width,
                        colorbar=false,
                    )
                else
                    error("Unsupported projection.plot_style: $(config.projection.plot_style)")
                end
            end
        end
    end

    for xpos in result.tick_positions
        vline!(plt, [xpos]; color=:gray, linewidth=0.8, linestyle=:dash)
    end
    xlabel!(plt, "")
    ylabel!(plt, config.band_plot.ylabel)
    xticks!(plt, result.tick_positions, result.tick_labels)
    xlims!(plt, extrema(result.distances))
    ylims!(plt, config.plot.energy_range)
    if save
        mkpath(dirname(config.projection.projected_bands_plot))
        savefig(plt, config.projection.projected_bands_plot)
    end
    return plt
end

function plot_dos(result::DosResult, config::RunConfig; save::Bool=true)
    plt = plot(
        size=config.plot.size,
        legend=false,
        dpi=200,
        framestyle=:box;
        _plot_grid_kwargs()...,
        _font_kwargs(config)...,
        _subplot_margins()...,
    )

    if !isnothing(result.dos_down)
        # result.dos makes a confusion with the non-spin polarized one.
        plot!(plt, result.energies, result.dos; color=config.spin.colors[1], linewidth=config.dos_plot.linewidth)
        plot!(plt, result.energies, -result.dos_down; color=config.spin.colors[2], linewidth=config.dos_plot.linewidth)
        hline!(plt, [0.0]; color=:gray, linewidth=0.6, linestyle=:dash)
    else
        plot!(plt, result.energies, result.dos; color=config.dos_plot.color, linewidth=config.dos_plot.linewidth)
    end

    _plot_center_vline!(plt, result.center_of_mass, config)
    xlabel!(plt, config.dos_plot.xlabel)
    ylabel!(plt, config.dos_plot.ylabel)
    xlims!(plt, config.plot.energy_range)
    if save
        mkpath(dirname(config.output.dos_plot))
        savefig(plt, config.output.dos_plot)
    end
    return plt
end

function plot_projected_dos(result::DosResult, config::RunConfig; save::Bool=true)
    projected = result.projected
    isnothing(projected) && error("DosResult has no projected DOS data")

    plt = plot(
        size=config.plot.size,
        legend=false,
        dpi=200,
        framestyle=:box;
        _plot_grid_kwargs()...,
        _font_kwargs(config)...,
        _subplot_margins()...,
    )

    for (ig, label) in enumerate(projected.labels)
        plot!(
            plt,
            result.energies,
            projected.pdos[:, ig];
            color=projected.colors[ig],
            linewidth=config.dos_plot.linewidth,
            label=label,
        )
        _plot_center_vline!(plt, projected.centers_of_mass[ig], config; color=projected.colors[ig])
    end

    xlabel!(plt, config.dos_plot.xlabel)
    ylabel!(plt, "PDOS (states/eV)")
    xlims!(plt, config.plot.energy_range)
    if save
        mkpath(dirname(config.projection.pdos_plot))
        savefig(plt, config.projection.pdos_plot)
    end
    return plt
end

function plot_combined(bands::BandResult, dos::DosResult, config::RunConfig)
    p1 = plot_bands(bands, config; save=false)

    p2 = plot(
        legend=false,
        dpi=200,
        framestyle=:box,
        xticks=_hidden_ticks(),
        yticks=[];
        _plot_grid_kwargs()...,
        _font_kwargs(config)...,
        _subplot_margins()...,
    )

    if !isnothing(dos.dos_down)
        plot!(p2, dos.dos, dos.energies; color=config.spin.colors[1], linewidth=config.dos_plot.linewidth)
        plot!(p2, -dos.dos_down, dos.energies; color=config.spin.colors[2], linewidth=config.dos_plot.linewidth)
        vline!(p2, [0.0]; color=:gray, linewidth=0.6, linestyle=:dash)
    else
        plot!(p2, dos.dos, dos.energies; color=config.dos_plot.color, linewidth=config.dos_plot.linewidth)
    end

    _plot_center_hline!(p2, dos.center_of_mass, config)
    xlabel!(p2, config.combined_plot.dos_xlabel)
    ylabel!(p2, config.combined_plot.dos_ylabel)
    ylims!(p2, config.plot.energy_range)
    plt = plot(
        p1,
        p2;
        layout=grid(1, 2; widths=_combined_widths(config)),
        size=(round(Int, config.plot.size[1] * (1 + config.combined_plot.dos_width_ratio)), config.plot.size[2]),
        _plot_grid_kwargs()...,
        _font_kwargs(config)...,
    )
    mkpath(dirname(config.output.combined_plot))
    savefig(plt, config.output.combined_plot)
    return plt
end

function plot_projected_combined(bands::BandResult, dos::DosResult, config::RunConfig)
    p1 = plot_projected_bands(bands, config; save=false)
    projected = dos.projected
    isnothing(projected) && error("DosResult has no projected DOS data")

    p2 = plot(
        legend=false,
        dpi=200,
        framestyle=:box,
        xticks=_hidden_ticks(),
        yticks=[];
        _plot_grid_kwargs()...,
        _font_kwargs(config)...,
        _subplot_margins()...,
    )

    for (ig, label) in enumerate(projected.labels)
        plot!(
            p2,
            projected.pdos[:, ig],
            dos.energies;
            color=projected.colors[ig],
            linewidth=config.dos_plot.linewidth,
            label=label,
        )
        _plot_center_hline!(p2, projected.centers_of_mass[ig], config; color=projected.colors[ig])
    end

    xlabel!(p2, config.combined_plot.dos_xlabel)
    ylabel!(p2, config.combined_plot.dos_ylabel)
    ylims!(p2, config.plot.energy_range)
    plt = plot(
        p1,
        p2;
        layout=grid(1, 2; widths=_combined_widths(config)),
        size=(round(Int, config.plot.size[1] * (1 + config.combined_plot.dos_width_ratio)), config.plot.size[2]),
        _plot_grid_kwargs()...,
        _font_kwargs(config)...,
    )
    mkpath(dirname(config.projection.projected_combined_plot))
    savefig(plt, config.projection.projected_combined_plot)
    return plt
end

end
