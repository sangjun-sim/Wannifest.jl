module Surface

using ..Execution
using ..Mesh: PlaneMesh, grid_index
using ..Model: EnergySurfaceResult, RunConfig
using ..WannierEigensystem

export compute_energy_surface

function compute_energy_surface(
    hr,
    mesh::PlaneMesh,
    config::RunConfig;
    wsvec=nothing,
)::EnergySurfaceResult
    nb = length(config.energy.bands)
    energies = Array{Float64, 3}(undef, mesh.ny, mesh.nx, nb)

    Execution.with_blas_threads(1) do
        Base.Threads.@threads for ik in eachindex(mesh.kpoints)
            evals, _ = WannierEigensystem.solve_kpoint(
                hr,
                mesh.kpoints[ik];
                wsvec=wsvec,
                hermiticity_tol=config.hermiticity_tol,
            )

            idx = grid_index(mesh, ik)
            for (out_band, band) in enumerate(config.energy.bands)
                1 <= band <= length(evals) ||
                    error("band index $band is outside 1:$(length(evals))")
                energies[idx.iy, idx.ix, out_band] = evals[band] - config.energy.shift
            end
        end
    end

    return EnergySurfaceResult(
        mesh.x_axis,
        mesh.y_axis,
        mesh.kpoints,
        copy(config.energy.bands),
        energies,
    )
end

end
