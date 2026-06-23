module Service

using Printf

using ..Mesh
using ..Model: ContourRunResult, RunConfig
using ..Output
using ..PlotService
using ..Surface
using ..WannierHrIO
using ..WannierWsvecGenerate
using ..WannierWsvecIO

export run, print_summary

_validate_wsvec(hr, wsvec) = WannierWsvecGenerate.assert_wsvec_usable(hr, wsvec)

function _check_hr_hermiticity(hr; atol::Float64=1e-10)
    maxdiff = WannierHrIO.pair_hermiticity_error(hr)
    if maxdiff > atol
        @warn "hr.dat deviates from Hermitian-pair symmetry by $maxdiff (threshold: $atol)"
        return false
    end
    return true
end

function run(config::RunConfig; make_plot::Bool=true)::ContourRunResult
    hr = WannierHrIO.read_hr(config.files.hr_path; spin_layout=config.spin_layout)
    hermiticity_ok = _check_hr_hermiticity(hr; atol=config.hermiticity_tol)
    wsvec = isnothing(config.files.wsvec_path) ? nothing : WannierWsvecIO.read_wsvec(
        config.files.wsvec_path;
        num_wann=hr.num_wann,
        spin_layout=config.spin_layout,
    )
    isnothing(wsvec) || _validate_wsvec(hr, wsvec)

    mesh = Mesh.generate_plane_mesh(config.plane)
    surface = Surface.compute_energy_surface(hr, mesh, config; wsvec=wsvec)

    data_path = Output.default_data_path(config)
    Output.write_surface_data(data_path, surface, config)
    plot_handles = PlotService.maybe_plot(config, surface, make_plot)

    return ContourRunResult(config, surface, data_path, plot_handles, hermiticity_ok)
end

function print_summary(result::ContourRunResult; make_plot::Bool=true, io::IO=stdout)
    println(io, "Contour run complete.")
    println(io, "Hermiticity check: ", result.hermiticity_ok ? "ok" : "failed")
    println(io, "Data: ", result.data_path)
    println(io, "Bands: ", join(result.surface.bands, ", "))
    @printf(io, "Mesh: %d x %d\n", length(result.surface.x_axis), length(result.surface.y_axis))
    if make_plot
        println(io, "Plots: ", length(result.plot_handles))
    end
    return nothing
end

end
