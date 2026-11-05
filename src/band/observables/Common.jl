module ObservablesCommon

using ..LocalAxisRotation
using ..Model: RunConfig
using ..Projection
using ..Win90Basis

export basis_entries_from_projection, expectation_values, warn_degenerate

function basis_entries_from_projection(config::RunConfig, num_wann::Integer, observable::AbstractString)
    projection = config.projection
    projection.enabled ||
        error("band.$observable requires band.projection metadata with mode=\"win_groups\"")
    projection.mode == :win_groups ||
        error("band.$observable cannot use band.projection mode=\"$(projection.mode)\"; use win_groups")
    isnothing(projection.win_path) && error("band.$observable requires band.projection.win")
    basis = Win90Basis.read_win_basis(projection.win_path; spin_layout=config.spin.layout)
    basis.num_wann == num_wann ||
        error("$(uppercase(observable)) basis has num_wann=$(basis.num_wann), but hr has num_wann=$num_wann")
    entries = LocalAxisRotation.basis_entries_from_win(basis)
    return Projection._canonicalize_basis_entries(entries, num_wann, config.spin.layout)
end

function expectation_values(expectations::Function, context, eigenvectors)
    if isnothing(context.basis_transform)
        return expectations(context.operators, eigenvectors)
    end
    size(context.basis_transform, 2) == size(eigenvectors, 1) ||
        error("basis rotation has $(size(context.basis_transform, 2)) columns, but evecs has $(size(eigenvectors, 1)) rows")
    return expectations(context.operators, context.basis_transform * eigenvectors)
end

function warn_degenerate(
    observable::AbstractString,
    eigenvalues::AbstractVector{<:Real},
    degeneracy_tol::Float64,
    ik::Integer,
    counter::Base.Threads.Atomic{Int};
    max_warnings::Int=10,
)
    for ib in 1:(length(eigenvalues) - 1)
        delta = abs(Float64(eigenvalues[ib + 1]) - Float64(eigenvalues[ib]))
        delta < degeneracy_tol || continue
        old = Base.Threads.atomic_add!(counter, 1)
        old < max_warnings && @warn(
            "band.$observable per-band values are gauge-dependent for degenerate states " *
            "(ik=$ik, bands=$ib/$(ib + 1), ΔE=$delta < degeneracy_tol=$degeneracy_tol)."
        )
        return nothing
    end
    return nothing
end

end
