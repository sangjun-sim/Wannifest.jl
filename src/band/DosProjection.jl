function _projected_centers_of_mass(energies::Vector{Float64}, pdos::Matrix{Float64})
    centers = Vector{Float64}(undef, size(pdos, 2))
    for ig in axes(pdos, 2)
        centers[ig] = density_center_of_mass(energies, view(pdos, :, ig))
    end
    return centers
end

function _accumulate_dos!(
    buffer::AbstractVector{Float64},
    evals::AbstractVector{<:Real},
    energies::AbstractVector{<:Real},
    sigma::Real,
)
    for eps in evals
        ilo = searchsortedfirst(energies, eps - 6 * sigma)
        ihi = searchsortedlast(energies, eps + 6 * sigma)
        ilo <= ihi || continue
        for ie in ilo:ihi
            buffer[ie] += gaussian(energies[ie], eps, sigma)
        end
    end
    return buffer
end

function _sum_vector_buffers(buffers)
    total = zeros(Float64, length(first(buffers)))
    for buffer in buffers
        total .+= buffer
    end
    return total
end

function _sum_matrix_buffers(buffers)
    total = zeros(Float64, size(first(buffers))...)
    for buffer in buffers
        total .+= buffer
    end
    return total
end

function _explicit_energy_grid(config::RunConfig)
    (isnothing(config.dos.emin) || isnothing(config.dos.emax)) && return nothing
    config.dos.emin < config.dos.emax || error("DOS energy window must have emin < emax")
    return collect(range(config.dos.emin, config.dos.emax; length=config.dos.npts))
end

function _dos_result_from_density(
    energies::Vector{Float64},
    dos::Vector{Float64},
    dos_down::Union{Nothing, Vector{Float64}},
    nbands::Int,
    window_is_auto::Bool;
    projected::Union{Nothing, ProjectedDosResult}=nothing,
)
    total_dos = isnothing(dos_down) ? dos : dos .+ dos_down
    integral = _check_dos_integral(total_dos, energies, nbands, window_is_auto)
    center_of_mass = density_center_of_mass(energies, total_dos)
    return DosResult(energies, dos, dos_down, integral, center_of_mass, nbands, window_is_auto, projected)
end

function _run_projected_dos_explicit(
    hr::HrData,
    kmesh::Vector{Vector{Float64}},
    wsvec::Union{Nothing, WsvecTable},
    config::RunConfig,
    spec::ProjectionSpec,
    basis_transform::Union{Nothing, AbstractMatrix{<:Complex}},
    energies::Vector{Float64},
)
    nthreads = Base.Threads.maxthreadid()
    ngroups = length(spec.groups)
    use_wsvec = !isnothing(wsvec)
    projected_buffers = [zeros(Float64, length(energies), ngroups) for _ in 1:nthreads]
    dos_buffers = [zeros(Float64, length(energies)) for _ in 1:nthreads]
    dos_down_buffers = config.spin.enabled ? [zeros(Float64, length(energies)) for _ in 1:nthreads] : nothing
    up_indices, down_indices = config.spin.enabled ?
        SpinLayout.band_indices(hr.num_wann, config.spin.layout) :
        (Int[], Int[])

    Execution.with_blas_threads(1) do
        Base.Threads.@threads for ik in eachindex(kmesh)
            evals, evecs = WannierEigensystem.solve_kpoint(
                hr,
                kmesh[ik];
                wsvec=use_wsvec ? wsvec : nothing,
                hermiticity_tol=config.hermiticity_tol,
            )
            shifted = evals .- config.energy.shift
            tid = Base.Threads.threadid()
            if config.spin.enabled
                _accumulate_dos!(dos_buffers[tid], view(shifted, up_indices), energies, config.dos.sigma)
                _accumulate_dos!(dos_down_buffers[tid], view(shifted, down_indices), energies, config.dos.sigma)
            else
                _accumulate_dos!(dos_buffers[tid], shifted, energies, config.dos.sigma)
            end
            weights = Projection.projection_weights(spec, evecs; basis_transform=basis_transform)
            accumulate_projected_dos!(
                projected_buffers[tid],
                shifted,
                weights,
                energies,
                config.dos.sigma,
            )
        end
    end

    nk = length(kmesh)
    dos = _sum_vector_buffers(dos_buffers) ./ nk
    dos_down = isnothing(dos_down_buffers) ? nothing : _sum_vector_buffers(dos_down_buffers) ./ nk
    pdos = _sum_matrix_buffers(projected_buffers) ./ nk
    projected = Projection.projected_dos_result(
        spec,
        pdos;
        centers_of_mass=_projected_centers_of_mass(energies, pdos),
    )
    return _dos_result_from_density(energies, dos, dos_down, hr.num_wann, false; projected=projected)
end

function _projected_dos_result(
    hr::HrData,
    kmesh::Vector{Vector{Float64}},
    wsvec::Union{Nothing, WsvecTable},
    config::RunConfig,
    spec::ProjectionSpec,
    energies::Vector{Float64},
    basis_transform::Union{Nothing, AbstractMatrix{<:Complex}},
)
    ngroups = length(spec.groups)
    buffers = [zeros(Float64, length(energies), ngroups) for _ in 1:Base.Threads.maxthreadid()]
    use_wsvec = !isnothing(wsvec)
    Execution.with_blas_threads(1) do
        Base.Threads.@threads for ik in eachindex(kmesh)
            evals, evecs = WannierEigensystem.solve_kpoint(
                hr,
                kmesh[ik];
                wsvec=use_wsvec ? wsvec : nothing,
                hermiticity_tol=config.hermiticity_tol,
            )
            weights = Projection.projection_weights(spec, evecs; basis_transform=basis_transform)
            accumulate_projected_dos!(
                buffers[Base.Threads.threadid()],
                evals .- config.energy.shift,
                weights,
                energies,
                config.dos.sigma,
            )
        end
    end

    pdos = _sum_matrix_buffers(buffers) ./ length(kmesh)
    centers_of_mass = _projected_centers_of_mass(energies, pdos)
    return Projection.projected_dos_result(spec, pdos; centers_of_mass=centers_of_mass)
end
