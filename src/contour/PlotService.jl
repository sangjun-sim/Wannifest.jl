module PlotService

const CONTOUR_DIR = @__DIR__

using ..Model: EnergySurfaceResult, RunConfig

export maybe_plot, plotcontour_module

const _plotcontour_module_ref = Ref{Union{Nothing, Module}}(nothing)

function _plotcontour_module()
    cached = _plotcontour_module_ref[]
    !isnothing(cached) && return cached

    parent = parentmodule(@__MODULE__)
    if !Base.invokelatest(isdefined, parent, :PlotContour)
        Base.include(parent, joinpath(CONTOUR_DIR, "PlotContour.jl"))
    end
    mod = Base.invokelatest(getfield, parent, :PlotContour)
    _plotcontour_module_ref[] = mod
    return mod
end

plotcontour_module() = _plotcontour_module()

function maybe_plot(config::RunConfig, result::EnergySurfaceResult, make_plot::Bool)
    make_plot || return Any[]
    plot_module = _plotcontour_module()
    handles = Any[]

    if config.plot.mode in (:surface, :both)
        fn = Base.invokelatest(getfield, plot_module, :plot_surface)
        push!(handles, Base.invokelatest(fn, result, config; save=true))
    end
    if config.plot.mode in (:contour, :both)
        fn = Base.invokelatest(getfield, plot_module, :plot_contour)
        push!(handles, Base.invokelatest(fn, result, config; save=true))
    end
    if config.plot.mode == :heatmap
        fn = Base.invokelatest(getfield, plot_module, :plot_heatmap)
        push!(handles, Base.invokelatest(fn, result, config; save=true))
    end

    if config.plot.interactive
        fn = Base.invokelatest(getfield, plot_module, :open_interactive)
        Base.invokelatest(fn, result, config)
    end
    return handles
end

end
