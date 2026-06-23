module PlotService

const BAND_DIR = @__DIR__

using ..Model: BandResult, DosResult, RunConfig

export maybe_plot, plotbands_module, observables_plot_module

const _plotbands_module_ref = Ref{Union{Nothing, Module}}(nothing)
const _observables_plot_module_ref = Ref{Union{Nothing, Module}}(nothing)

function _plotbands_module()
    cached = _plotbands_module_ref[]
    !isnothing(cached) && return cached

    parent = parentmodule(@__MODULE__)
    if !Base.invokelatest(isdefined, parent, :PlotBands)
        Base.include(parent, joinpath(BAND_DIR, "PlotBands.jl"))
    end
    mod = Base.invokelatest(getfield, parent, :PlotBands)
    _plotbands_module_ref[] = mod
    return mod
end

plotbands_module() = _plotbands_module()

function _observables_plot_module()
    cached = _observables_plot_module_ref[]
    !isnothing(cached) && return cached

    parent = parentmodule(@__MODULE__)
    if !Base.invokelatest(isdefined, parent, :ObservablesPlot)
        Base.include(parent, joinpath(BAND_DIR, "observables", "Plot.jl"))
    end
    mod = Base.invokelatest(getfield, parent, :ObservablesPlot)
    _observables_plot_module_ref[] = mod
    return mod
end

observables_plot_module() = _observables_plot_module()

_wants_plot(config::RunConfig, target::Symbol) = target in config.plot.targets

function maybe_plot(
    config::RunConfig,
    band_result::Union{Nothing, BandResult},
    dos_result::Union{Nothing, DosResult},
    make_plot::Bool,
)
    make_plot || return nothing
    isempty(config.plot.targets) && return nothing
    interactive_plt = nothing

    if !isnothing(band_result) && _wants_plot(config, :band)
        plot_module = _plotbands_module()
        plot_bands_fn = Base.invokelatest(getfield, plot_module, :plot_bands)
        plt = Base.invokelatest(plot_bands_fn, band_result, config)
        isnothing(interactive_plt) && (interactive_plt = plt)
    end

    if !isnothing(band_result) && config.projection.enabled &&
            !isnothing(band_result.projection) && _wants_plot(config, :fatband)
        plot_module = _plotbands_module()
        plot_projected_bands_fn = Base.invokelatest(getfield, plot_module, :plot_projected_bands)
        plt = Base.invokelatest(plot_projected_bands_fn, band_result, config)
        isnothing(interactive_plt) && (interactive_plt = plt)
    end

    if !isnothing(band_result) && config.oam.enabled &&
            !isnothing(band_result.oam) && _wants_plot(config, :oam)
        oam_plot_module = _observables_plot_module()
        plot_oam_bands_fn = Base.invokelatest(getfield, oam_plot_module, :plot_oam_bands)
        for component in config.oam.plot_components
            plt = Base.invokelatest(plot_oam_bands_fn, band_result, config; component=component)
            isnothing(interactive_plt) && (interactive_plt = plt)
        end
    end

    if !isnothing(band_result) && config.sam.enabled &&
            !isnothing(band_result.sam) && _wants_plot(config, :sam)
        sam_plot_module = _observables_plot_module()
        plot_sam_bands_fn = Base.invokelatest(getfield, sam_plot_module, :plot_sam_bands)
        for component in config.sam.plot_components
            plt = Base.invokelatest(plot_sam_bands_fn, band_result, config; component=component)
            isnothing(interactive_plt) && (interactive_plt = plt)
        end
    end

    if !isnothing(dos_result) && _wants_plot(config, :dos)
        plot_module = _plotbands_module()
        plot_dos_fn = Base.invokelatest(getfield, plot_module, :plot_dos)
        plt = Base.invokelatest(plot_dos_fn, dos_result, config)
        isnothing(interactive_plt) && (interactive_plt = plt)
    end

    if !isnothing(dos_result) && config.projection.enabled &&
            !isnothing(dos_result.projected) && _wants_plot(config, :pdos)
        plot_module = _plotbands_module()
        plot_projected_dos_fn = Base.invokelatest(getfield, plot_module, :plot_projected_dos)
        plt = Base.invokelatest(plot_projected_dos_fn, dos_result, config)
        isnothing(interactive_plt) && (interactive_plt = plt)
    end

    if !isnothing(band_result) && !isnothing(dos_result) && _wants_plot(config, :combined)
        plot_module = _plotbands_module()
        plot_combined_fn = Base.invokelatest(getfield, plot_module, :plot_combined)
        plt = Base.invokelatest(plot_combined_fn, band_result, dos_result, config)
        isnothing(interactive_plt) && (interactive_plt = plt)
    end

    if !isnothing(band_result) && !isnothing(dos_result) && config.projection.enabled &&
            !isnothing(band_result.projection) && !isnothing(dos_result.projected) &&
            _wants_plot(config, :fatband_pdos)
        plot_module = _plotbands_module()
        plot_projected_combined_fn = Base.invokelatest(getfield, plot_module, :plot_projected_combined)
        plt = Base.invokelatest(plot_projected_combined_fn, band_result, dos_result, config)
        isnothing(interactive_plt) && (interactive_plt = plt)
    end

    if config.plot.interactive && !isnothing(interactive_plt)
        Base.invokelatest(display, interactive_plt)
    end

    return nothing
end

end
