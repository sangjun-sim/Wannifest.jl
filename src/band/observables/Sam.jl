module ObservablesSam

using ..AtomicSpin
using ..Model: RunConfig
using ..ObservablesCommon

export SamContext, build_sam_context, sam_expectation_values, warn_degenerate_sam

struct SamContext
    operators::AtomicSpin.SpinOperators
    basis_transform::Union{Nothing, Matrix{ComplexF64}}
end

function build_sam_context(
    config::RunConfig,
    num_wann::Integer;
    projection_basis_transform::Union{Nothing, AbstractMatrix{<:Complex}}=nothing,
)::Union{Nothing, SamContext}
    config.sam.enabled || return nothing
    entries = basis_entries_from_projection(config, num_wann, "sam")
    operators = AtomicSpin.build_spin_operators(num_wann, entries)
    transform = isnothing(projection_basis_transform) ? nothing : Matrix{ComplexF64}(projection_basis_transform)
    return SamContext(operators, transform)
end

sam_expectation_values(context::SamContext, eigenvectors) =
    expectation_values(AtomicSpin.spin_expectations, context, eigenvectors)

warn_degenerate_sam(args...; kwargs...) = warn_degenerate("sam", args...; kwargs...)

end
