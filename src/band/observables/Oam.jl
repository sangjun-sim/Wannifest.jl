module ObservablesOam

using ..AtomicOam
using ..Model: RunConfig
using ..ObservablesCommon

export OamContext, build_oam_context, oam_expectation_values, warn_degenerate_oam

struct OamContext
    operators::AtomicOam.OamOperators
    basis_transform::Union{Nothing, Matrix{ComplexF64}}
end

function build_oam_context(
    config::RunConfig,
    num_wann::Integer;
    projection_basis_transform::Union{Nothing, AbstractMatrix{<:Complex}}=nothing,
)::Union{Nothing, OamContext}
    config.oam.enabled || return nothing
    entries = basis_entries_from_projection(config, num_wann, "oam")
    operators = AtomicOam.build_l_operators(num_wann, entries, config.oam.selections)
    transform = isnothing(projection_basis_transform) ? nothing : Matrix{ComplexF64}(projection_basis_transform)
    return OamContext(operators, transform)
end

oam_expectation_values(context::OamContext, eigenvectors) =
    expectation_values(AtomicOam.oam_expectations, context, eigenvectors)

warn_degenerate_oam(args...; kwargs...) = warn_degenerate("oam", args...; kwargs...)

end
