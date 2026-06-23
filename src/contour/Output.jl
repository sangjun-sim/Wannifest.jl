module Output

using Printf

using ..Model: EnergySurfaceResult, RunConfig

export default_data_path, default_plot_path, default_interactive_path, write_surface_data

function default_data_path(config::RunConfig)
    return joinpath(config.output.output_dir, "data", "contour_energy_surface.dat")
end

function default_plot_path(config::RunConfig, kind::Symbol)
    filename = "contour_$(kind).png"
    return joinpath(config.output.output_dir, "plots", filename)
end

function default_interactive_path(config::RunConfig, kind::Symbol)
    filename = "contour_$(kind).html"
    return joinpath(config.output.output_dir, "plots", filename)
end

function write_surface_data(path::AbstractString, result::EnergySurfaceResult, config::RunConfig)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# columns = ix iy kx ky kz band energy")
        println(io, "# energy_shift = ", config.energy.shift)
        println(io, "# hr = ", abspath(config.files.hr_path))
        println(io, "# structure = ", isempty(config.files.structure_path) ? "" : abspath(config.files.structure_path))
        println(io, "# wsvec = ", isnothing(config.files.wsvec_path) ? "" : abspath(config.files.wsvec_path))
        println(io, "# spin_layout = ", config.spin_layout)
        for (ib_out, band) in enumerate(result.bands)
            for iy in eachindex(result.y_axis), ix in eachindex(result.x_axis)
                ik = (iy - 1) * length(result.x_axis) + ix
                k = result.kpoints_frac[ik]
                @printf(
                    io,
                    "%d %d %.10f %.10f %.10f %d %.10f\n",
                    ix,
                    iy,
                    k[1],
                    k[2],
                    k[3],
                    band,
                    result.energies[iy, ix, ib_out],
                )
            end
        end
    end
    return path
end

end
