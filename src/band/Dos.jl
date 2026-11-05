module Dos

using ..Execution
using ..Model: DosResult, HrData, RunConfig, WsvecTable
using ..Model: _density_center_of_mass, _trapezoidal_integral
using ..OrbitalProjection: ProjectionSpec
using ..Projection
using ..ProjectionModel: ProjectedDosResult
using ..SpinLayout
using ..WannierEigensystem

export generate_kmesh, compute_dos, run_dos, gaussian, accumulate_projected_dos!
export trapezoidal_integral, density_center_of_mass

trapezoidal_integral(values::AbstractVector{<:Real}, coordinates::AbstractVector{<:Real}) =
    _trapezoidal_integral(values, coordinates)

density_center_of_mass(energies::AbstractVector{<:Real}, density::AbstractVector{<:Real}) =
    _density_center_of_mass(energies, density)

function generate_kmesh(n1::Int, n2::Int, n3::Int; shift::NTuple{3, Float64}=(0.0, 0.0, 0.0))
    n1 > 0 && n2 > 0 && n3 > 0 || error("k-mesh dimensions must be positive")
    all(x -> 0.0 <= x < 1.0, shift) || error("k-mesh shift entries must satisfy 0.0 <= shift < 1.0")
    kpoints = Vector{Float64}[]
    for i3 in 0:(n3 - 1), i2 in 0:(n2 - 1), i1 in 0:(n1 - 1)
        push!(kpoints, [
            mod((i1 + shift[1]) / n1, 1.0),
            mod((i2 + shift[2]) / n2, 1.0),
            mod((i3 + shift[3]) / n3, 1.0),
        ])
    end
    return kpoints
end

function _gaussian_dos(eigenvalues::AbstractMatrix{Float64}, energies::Vector{Float64}, sigma::Float64)::Vector{Float64}
    nk, nbands = size(eigenvalues)
    npts = length(energies)
    dos = zeros(Float64, npts)
    prefactor = (1.0 / nk) / (sigma * sqrt(2π))
    local_buffers = [zeros(Float64, npts) for _ in 1:Base.Threads.maxthreadid()]

    Base.Threads.@threads for ik in 1:nk
        tid = Base.Threads.threadid()
        buf = local_buffers[tid]
        for ib in 1:nbands
            eps = eigenvalues[ik, ib]
            ilo = searchsortedfirst(energies, eps - 6 * sigma)
            ihi = searchsortedlast(energies, eps + 6 * sigma)
            ilo <= ihi || continue
            for ie in ilo:ihi
                x = (energies[ie] - eps) / sigma
                buf[ie] += prefactor * exp(-0.5 * x * x)
            end
        end
    end

    for buf in local_buffers
        dos .+= buf
    end
    return dos
end

function gaussian(e::Real, e0::Real, sigma::Real)
    x = (Float64(e) - Float64(e0)) / Float64(sigma)
    return exp(-0.5 * x * x) / (Float64(sigma) * sqrt(2π))
end

function accumulate_projected_dos!(
    buffer::AbstractMatrix{Float64},
    evals::AbstractVector{<:Real},
    weights::AbstractMatrix{Float64},
    energies::AbstractVector{<:Real},
    sigma::Real,
)
    size(weights, 1) == length(evals) ||
        error("projected DOS weights have $(size(weights, 1)) bands, expected $(length(evals))")
    size(buffer, 1) == length(energies) ||
        error("projected DOS buffer has $(size(buffer, 1)) energies, expected $(length(energies))")
    size(buffer, 2) == size(weights, 2) ||
        error("projected DOS buffer has $(size(buffer, 2)) groups, expected $(size(weights, 2))")

    for (ib, eps) in enumerate(evals)
        ilo = searchsortedfirst(energies, eps - 6 * sigma)
        ihi = searchsortedlast(energies, eps + 6 * sigma)
        ilo <= ihi || continue
        for ig in axes(weights, 2)
            weight = weights[ib, ig]
            weight == 0.0 && continue
            for ie in ilo:ihi
                buffer[ie, ig] += weight * gaussian(energies[ie], eps, sigma)
            end
        end
    end
    return buffer
end

function _check_dos_integral(dos::Vector{Float64}, energies::Vector{Float64}, nbands::Int, window_is_auto::Bool)
    integral = trapezoidal_integral(dos, energies)
    deviation = abs(integral - nbands) / nbands
    if window_is_auto && deviation > 0.01
        @warn "DOS integral = $integral, expected $nbands (deviation $(round(deviation * 100; digits=2))%). Consider increasing mesh density or adjusting sigma."
    elseif !window_is_auto && deviation > 0.01
        @info "DOS integral over requested energy window is $integral; full-spectrum normalization would be $nbands."
    end
    return integral
end

function compute_dos(
    all_eigenvalues::Matrix{Float64},
    config::RunConfig;
    projected::Union{Nothing, ProjectedDosResult}=nothing,
)
    nbands = size(all_eigenvalues, 2)
    sigma = config.dos.sigma
    window_is_auto = isnothing(config.dos.emin) && isnothing(config.dos.emax)
    emin = isnothing(config.dos.emin) ? minimum(all_eigenvalues) - 5 * sigma : config.dos.emin
    emax = isnothing(config.dos.emax) ? maximum(all_eigenvalues) + 5 * sigma : config.dos.emax
    emin < emax || error("DOS energy window must have emin < emax")
    npts = config.dos.npts
    energies = collect(range(emin, emax; length=npts))

    if config.spin.enabled
        up_indices, down_indices = SpinLayout.band_indices(nbands, config.spin.layout)
        dos_up = _gaussian_dos(view(all_eigenvalues, :, up_indices), energies, sigma)
        dos_dn = _gaussian_dos(view(all_eigenvalues, :, down_indices), energies, sigma)
        total_dos = dos_up .+ dos_dn
        integral = _check_dos_integral(total_dos, energies, nbands, window_is_auto)
        center_of_mass = density_center_of_mass(energies, total_dos)
        return DosResult(energies, dos_up, dos_dn, integral, center_of_mass, nbands, window_is_auto, projected)
    end

    dos = _gaussian_dos(all_eigenvalues, energies, sigma)
    integral = _check_dos_integral(dos, energies, nbands, window_is_auto)
    center_of_mass = density_center_of_mass(energies, dos)
    return DosResult(energies, dos, nothing, integral, center_of_mass, nbands, window_is_auto, projected)
end

include(joinpath(@__DIR__, "DosProjection.jl"))

function run_dos(
    hr::HrData,
    config::RunConfig;
    wsvec::Union{Nothing, WsvecTable}=nothing,
    projection_spec::Union{Nothing, ProjectionSpec}=nothing,
    basis_transform::Union{Nothing, AbstractMatrix{<:Complex}}=nothing,
)
    kmesh = generate_kmesh(config.dos.mesh...; shift=config.dos.shift)
    if !isnothing(projection_spec)
        explicit_energies = _explicit_energy_grid(config)
        if !isnothing(explicit_energies)
            return _run_projected_dos_explicit(
                hr,
                kmesh,
                wsvec,
                config,
                projection_spec,
                basis_transform,
                explicit_energies,
            )
        end
    end

    all_evals = Matrix{Float64}(undef, length(kmesh), hr.num_wann)
    use_wsvec = !isnothing(wsvec)
    Execution.with_blas_threads(1) do
        Base.Threads.@threads for ik in eachindex(kmesh)
            evals = WannierEigensystem.solve_kpoint_values(
                hr,
                kmesh[ik];
                wsvec=use_wsvec ? wsvec : nothing,
                hermiticity_tol=config.hermiticity_tol,
            )
            all_evals[ik, :] .= evals .- config.energy.shift
        end
    end

    dos_result = compute_dos(all_evals, config)
    isnothing(projection_spec) && return dos_result

    projected = _projected_dos_result(
        hr,
        kmesh,
        wsvec,
        config,
        projection_spec,
        dos_result.energies,
        basis_transform,
    )
    return DosResult(
        dos_result.energies,
        dos_result.dos,
        dos_result.dos_down,
        dos_result.integral,
        dos_result.center_of_mass,
        dos_result.num_bands,
        dos_result.window_is_auto,
        projected,
    )
end

end
