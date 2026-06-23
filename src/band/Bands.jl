module Bands

using ..Execution
using ..KPath: KPathResult
using ..Model: BandResult, HrData, RunConfig, WsvecTable
using ..ObservablesOam
using ..ObservablesSam
using ..OrbitalProjection: ProjectionSpec
using ..Projection
using ..WannierEigensystem

export compute_bands

function compute_bands(
    hr::HrData,
    kpath::KPathResult,
    config::RunConfig;
    wsvec::Union{Nothing, WsvecTable}=nothing,
    projection_spec::Union{Nothing, ProjectionSpec}=nothing,
    basis_transform::Union{Nothing, AbstractMatrix{<:Complex}}=nothing,
    oam_context::Union{Nothing, ObservablesOam.OamContext}=nothing,
    sam_context::Union{Nothing, ObservablesSam.SamContext}=nothing,
)
    use_wsvec = !isnothing(wsvec)
    kpoints = kpath.kpoints
    band_eigenvalues = Matrix{Float64}(undef, length(kpoints), hr.num_wann)
    projection_weights = isnothing(projection_spec) ?
        nothing :
        Array{Float64, 3}(undef, length(kpoints), hr.num_wann, length(projection_spec.groups))
    oam_values = isnothing(oam_context) ? nothing : Array{Float64, 3}(undef, length(kpoints), hr.num_wann, 5)
    degeneracy_warnings = Base.Threads.Atomic{Int}(0)
    sam_values = isnothing(sam_context) ? nothing : Array{Float64, 3}(undef, length(kpoints), hr.num_wann, 5)
    sam_degeneracy_warnings = Base.Threads.Atomic{Int}(0)
    Execution.with_blas_threads(1) do
        Base.Threads.@threads for ik in eachindex(kpoints)
            evals, evecs = WannierEigensystem.solve_kpoint(
                hr,
                kpoints[ik];
                wsvec=use_wsvec ? wsvec : nothing,
                hermiticity_tol=config.hermiticity_tol,
            )
            band_eigenvalues[ik, :] .= evals
            if !isnothing(projection_weights) && !isnothing(projection_spec)
                projection_weights[ik, :, :] .= Projection.projection_weights(
                    projection_spec,
                    evecs;
                    basis_transform=basis_transform,
                )
            end
            if !isnothing(oam_values) && !isnothing(oam_context)
                oam_values[ik, :, :] .= ObservablesOam.oam_expectation_values(oam_context, evecs)
                ObservablesOam.warn_degenerate_oam(
                    evals,
                    config.oam.degeneracy_tol,
                    ik,
                    degeneracy_warnings,
                )
            end
            if !isnothing(sam_values) && !isnothing(sam_context)
                sam_values[ik, :, :] .= ObservablesSam.sam_expectation_values(sam_context, evecs)
                ObservablesSam.warn_degenerate_sam(
                    evals,
                    config.sam.degeneracy_tol,
                    ik,
                    sam_degeneracy_warnings,
                )
            end
        end
    end
    eigenvalues = band_eigenvalues .- config.energy.shift

    if !kpath.is_physical_distance && config.verbose
        @info "No lattice provided. Band plot x-axis uses reduced k-coordinate distance (not physical units)."
    end

    projection = if isnothing(projection_spec) || isnothing(projection_weights)
        nothing
    else
        Projection.band_projection_result(projection_spec, projection_weights)
    end

    return BandResult(
        kpoints,
        kpath.distances,
        eigenvalues,
        kpath.tick_positions,
        kpath.tick_labels,
        kpath.segment_ranges,
        kpath.is_physical_distance,
        projection,
        oam_values,
        sam_values,
    )
end

function compute_bands(
    hr::HrData,
    kpoints::Vector{Vector{Float64}},
    distances::Vector{Float64},
    tick_pos::Vector{Float64},
    tick_labels::Vector{String},
    segment_ranges::Vector{UnitRange{Int}},
    is_physical::Bool,
    config::RunConfig;
    wsvec::Union{Nothing, WsvecTable}=nothing,
    projection_spec::Union{Nothing, ProjectionSpec}=nothing,
    basis_transform::Union{Nothing, AbstractMatrix{<:Complex}}=nothing,
    oam_context::Union{Nothing, ObservablesOam.OamContext}=nothing,
    sam_context::Union{Nothing, ObservablesSam.SamContext}=nothing,
)
    kpath = KPathResult(kpoints, distances, tick_pos, tick_labels, segment_ranges, is_physical)
    return compute_bands(
        hr,
        kpath,
        config;
        wsvec=wsvec,
        projection_spec=projection_spec,
        basis_transform=basis_transform,
        oam_context=oam_context,
        sam_context=sam_context,
    )
end

end
