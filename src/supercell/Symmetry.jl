module Symmetry

using ..CellConventions

export SymmetrySummary
export build_spglib_cell
export build_dataset
export summarize_dataset

const SymmetrySummary = CellConventions.SymmetrySummary
build_spglib_cell(args...; kwargs...) = CellConventions.build_spglib_cell(args...; kwargs...)
build_dataset(args...; kwargs...) = CellConventions.build_dataset(args...; kwargs...)
summarize_dataset(args...; kwargs...) = CellConventions.summarize_dataset(args...; kwargs...)

end
