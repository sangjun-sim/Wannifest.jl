module PlotBands

using Plots

using ..Model: BandComparison, BandResult, DosResult, RunConfig
using ..Output
using ..SpinLayout

export plot_bands, plot_dos, plot_combined
export plot_projected_bands, plot_projected_dos, plot_projected_combined
export plot_compare

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

# An explicit plot.fermi_level always wins; otherwise the caller-supplied level
# (QE reference alignment) applies, and shifted energies default to zero.
_resolved_fermi_level(config::RunConfig, fermi_level::Union{Nothing, Real}) =
    Float64(something(config.plot.fermi_level, fermi_level, 0.0))

function _fermi_marker_kwargs()
    return (color=:gray, linewidth=0.9, linestyle=:dash, label=false)
end

function _plot_fermi_hline!(plt, config::RunConfig; fermi_level::Union{Nothing, Real}=nothing)
    config.plot.fermi_line || return plt
    level = _resolved_fermi_level(config, fermi_level)
    _center_visible(level, config) || return plt
    hline!(plt, [level]; _fermi_marker_kwargs()...)
    return plt
end

function _plot_fermi_vline!(plt, config::RunConfig; fermi_level::Union{Nothing, Real}=nothing)
    config.plot.fermi_line || return plt
    level = _resolved_fermi_level(config, fermi_level)
    _center_visible(level, config) || return plt
    vline!(plt, [level]; _fermi_marker_kwargs()...)
    return plt
end

# High-symmetry guide lines are drawn before any data series so the bands stay
# on top of them.
function _plot_tick_vlines!(plt, result::BandResult)
    for xpos in result.tick_positions
        vline!(plt, [xpos]; color=:gray, linewidth=0.8, linestyle=:dash, label=false)
    end
    return plt
end

function _band_canvas(config::RunConfig; legend=false)
    return plot(
        size=config.plot.size,
        legend=legend,
        dpi=200,
        framestyle=:box;
        _plot_grid_kwargs()...,
        _font_kwargs(config)...,
        _subplot_margins(config)...,
    )
end

function _finalize_band_axes!(
    plt,
    result::BandResult,
    config::RunConfig;
    fermi_level::Union{Nothing, Real}=nothing,
)
    _plot_fermi_hline!(plt, config; fermi_level=fermi_level)
    xlabel!(plt, "")
    ylabel!(plt, config.band_plot.ylabel)
    xticks!(plt, result.tick_positions, result.tick_labels)
    xlims!(plt, extrema(result.distances))
    ylims!(plt, config.plot.energy_range)
    return plt
end

function _save_plot(plt, path::Union{Nothing, AbstractString})
    isnothing(path) && return plt
    mkpath(dirname(path))
    savefig(plt, path)
    return plt
end

# Collinear runs (run.hr_down) get one figure per spin channel, reusing the
# `_up`/`_dn` suffixes that the band data files already carry.
_channel_plot_path(path::AbstractString, spin::Symbol) = Output.channel_path(path, spin)

function _plot_bands_channel(
    result::BandResult,
    config::RunConfig;
    color=nothing,
    label::Union{Nothing, String}=nothing,
    path::Union{Nothing, AbstractString}=nothing,
)
    show_legend = !isnothing(label) && config.band_plot.legend
    plt = _band_canvas(config; legend=show_legend ? :topright : false)

    _plot_tick_vlines!(plt, result)

    labeled = false
    nbands = size(result.eigenvalues, 2)
    for ib in _band_plot_order(config, nbands)
        for seg_range in result.segment_ranges
            plot!(
                plt,
                result.distances[seg_range],
                result.eigenvalues[seg_range, ib];
                color=isnothing(color) ? _band_color(config, ib, nbands) : color,
                linewidth=config.band_plot.linewidth,
                label=show_legend && !labeled ? label : false,
            )
            labeled = true
        end
    end

    _finalize_band_axes!(plt, result, config)
    return _save_plot(plt, path)
end

function plot_bands(
    result::BandResult,
    config::RunConfig;
    down::Union{Nothing, BandResult}=nothing,
    save::Bool=true,
)
    if isnothing(down)
        return _plot_bands_channel(
            result,
            config;
            path=save ? config.output.bands_plot : nothing,
        )
    end

    plt_up = _plot_bands_channel(
        result,
        config;
        color=config.spin.colors[1],
        label="up",
        path=save ? _channel_plot_path(config.output.bands_plot, :up) : nothing,
    )
    _plot_bands_channel(
        down,
        config;
        color=config.spin.colors[2],
        label="down",
        path=save ? _channel_plot_path(config.output.bands_plot, :down) : nothing,
    )
    return plt_up
end

function _plot_compare_channel(
    result::BandResult,
    comparison::BandComparison,
    config::RunConfig;
    wannier_color=nothing,
    qe_color=nothing,
    wannier_label::String="Wannier",
    qe_label::String="QE",
    fermi_level::Union{Nothing, Real}=nothing,
    path::Union{Nothing, AbstractString}=nothing,
)
    plt = _band_canvas(config; legend=:topright)

    _plot_tick_vlines!(plt, result)

    nbands = size(result.eigenvalues, 2)
    wannier_labeled = false
    qe_labeled = false
    for seg_range in result.segment_ranges
        for ib in _band_plot_order(config, nbands)
            plot!(
                plt,
                result.distances[seg_range],
                result.eigenvalues[seg_range, ib];
                color=isnothing(wannier_color) ? _band_color(config, ib, nbands) : wannier_color,
                linewidth=config.band_plot.linewidth,
                label=wannier_labeled ? false : wannier_label,
            )
            wannier_labeled = true
        end
        for iq in axes(comparison.qe_eigenvalues, 2)
            plot!(
                plt,
                result.distances[seg_range],
                comparison.qe_eigenvalues[seg_range, iq];
                color=isnothing(qe_color) ? config.reference.color : qe_color,
                linewidth=config.reference.linewidth,
                linestyle=:dash,
                label=qe_labeled ? false : qe_label,
            )
            qe_labeled = true
        end
    end

    _finalize_band_axes!(plt, result, config; fermi_level=fermi_level)
    return _save_plot(plt, path)
end

# Legend name of the DFT reference lines, after `[band.reference].source`.
_reference_label(config::RunConfig) = config.reference.source == :vasp ? "VASP" : "QE"

function plot_compare(
    result::BandResult,
    comparison::BandComparison,
    config::RunConfig;
    down::Union{Nothing, BandResult}=nothing,
    comparison_down::Union{Nothing, BandComparison}=nothing,
    fermi_level::Union{Nothing, Real}=nothing,
    save::Bool=true,
)
    reference_label = _reference_label(config)
    if isnothing(down)
        return _plot_compare_channel(
            result,
            comparison,
            config;
            qe_label=reference_label,
            fermi_level=fermi_level,
            path=save ? config.output.compare_plot : nothing,
        )
    end

    isnothing(comparison_down) && error("spin-down band result requires a spin-down reference comparison")
    plt_up = _plot_compare_channel(
        result,
        comparison,
        config;
        wannier_color=config.spin.colors[1],
        qe_color=config.reference.color,
        wannier_label="Wannier up",
        qe_label="$reference_label up",
        fermi_level=fermi_level,
        path=save ? _channel_plot_path(config.output.compare_plot, :up) : nothing,
    )
    _plot_compare_channel(
        down,
        comparison_down,
        config;
        wannier_color=config.spin.colors[2],
        qe_color=config.reference.color_down,
        wannier_label="Wannier down",
        qe_label="$reference_label down",
        fermi_level=fermi_level,
        path=save ? _channel_plot_path(config.output.compare_plot, :down) : nothing,
    )
    return plt_up
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

    plt = _band_canvas(config)

    _plot_tick_vlines!(plt, result)

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

    _plot_fermi_hline!(plt, config)
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
    plt = _band_canvas(config)

    if !isnothing(result.dos_down)
        # result.dos makes a confusion with the non-spin polarized one.
        plot!(plt, result.energies, result.dos; color=config.spin.colors[1], linewidth=config.dos_plot.linewidth)
        plot!(plt, result.energies, -result.dos_down; color=config.spin.colors[2], linewidth=config.dos_plot.linewidth)
        hline!(plt, [0.0]; color=:gray, linewidth=0.6, linestyle=:dash)
    else
        plot!(plt, result.energies, result.dos; color=config.dos_plot.color, linewidth=config.dos_plot.linewidth)
    end

    _plot_center_vline!(plt, result.center_of_mass, config)
    _plot_fermi_vline!(plt, config)
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

    plt = _band_canvas(config)

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

    _plot_fermi_vline!(plt, config)
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
        _subplot_margins(config)...,
    )

    if !isnothing(dos.dos_down)
        plot!(p2, dos.dos, dos.energies; color=config.spin.colors[1], linewidth=config.dos_plot.linewidth)
        plot!(p2, -dos.dos_down, dos.energies; color=config.spin.colors[2], linewidth=config.dos_plot.linewidth)
        vline!(p2, [0.0]; color=:gray, linewidth=0.6, linestyle=:dash)
    else
        plot!(p2, dos.dos, dos.energies; color=config.dos_plot.color, linewidth=config.dos_plot.linewidth)
    end

    _plot_center_hline!(p2, dos.center_of_mass, config)
    _plot_fermi_hline!(p2, config)
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
        _subplot_margins(config)...,
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

    _plot_fermi_hline!(p2, config)
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
