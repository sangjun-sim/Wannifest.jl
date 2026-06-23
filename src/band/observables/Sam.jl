module ObservablesSam

using ..AtomicSpin
using ..LocalAxisRotation
using ..Model: RunConfig
using ..Projection
using ..Win90Basis

export SamContext, build_sam_context, sam_expectation_values, warn_degenerate_sam

struct SamContext
    operators::AtomicSpin.SpinOperators
    basis_transform::Union{Nothing, Matrix{ComplexF64}}
end

function _basis_entries_from_projection(config::RunConfig)
    projection = config.projection
    projection.enabled ||
        error("band.sam requires band.projection metadata with mode=\"win_groups\"")
    if projection.mode == :win_groups
        isnothing(projection.win_path) && error("band.sam requires band.projection.win")
        basis = Win90Basis.read_win_basis(projection.win_path; spin_layout=config.spin.layout)
        return basis.num_wann, LocalAxisRotation.basis_entries_from_win(basis)
    end
    error("band.sam cannot use band.projection mode=\"$(projection.mode)\"; use win_groups")
end

function build_sam_context(
    config::RunConfig,
    num_wann::Integer;
    projection_basis_transform::Union{Nothing, AbstractMatrix{<:Complex}}=nothing,
)::Union{Nothing, SamContext}
    config.sam.enabled || return nothing
    basis_num_wann, entries = _basis_entries_from_projection(config)
    basis_num_wann == num_wann ||
        error("SAM basis has num_wann=$basis_num_wann, but hr has num_wann=$num_wann")
    canonical_entries = Projection._canonicalize_basis_entries(entries, num_wann, config.spin.layout)
    operators = AtomicSpin.build_spin_operators(num_wann, canonical_entries)
    transform = isnothing(projection_basis_transform) ? nothing : Matrix{ComplexF64}(projection_basis_transform)
    return SamContext(operators, transform)
end

function sam_expectation_values(context::SamContext, eigenvectors)
    if isnothing(context.basis_transform)
        return AtomicSpin.spin_expectations(context.operators, eigenvectors)
    end
    size(context.basis_transform, 2) == size(eigenvectors, 1) ||
        error("basis rotation has $(size(context.basis_transform, 2)) columns, but evecs has $(size(eigenvectors, 1)) rows")
    return AtomicSpin.spin_expectations(context.operators, context.basis_transform * eigenvectors)
end

function warn_degenerate_sam(
    eigenvalues::AbstractVector{<:Real},
    degeneracy_tol::Float64,
    ik::Integer,
    counter::Base.Threads.Atomic{Int};
    max_warnings::Int=10,
)
    nbands = length(eigenvalues)
    for ib in 1:(nbands - 1)
        delta = abs(Float64(eigenvalues[ib + 1]) - Float64(eigenvalues[ib]))
        delta < degeneracy_tol || continue
        old = Base.Threads.atomic_add!(counter, 1)
        old < max_warnings && @warn(
            "band.sam per-band values are gauge-dependent for degenerate states " *
            "(ik=$ik, bands=$ib/$(ib + 1), ΔE=$delta < degeneracy_tol=$degeneracy_tol)."
        )
        return nothing
    end
    return nothing
end

end
