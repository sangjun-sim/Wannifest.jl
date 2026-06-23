module Model

using ..WannierTypes: RKey

export FluxEndpoint, FluxRow, FluxTerm, RunFiles, GeometryConfig, ArrowStyle, PlotConfig, SpinConfig
export BasisConfig, DiagnosticVertex, PlaquetteDiagnostic, DiagnosticConfig, FluxConfig
export FluxBasisEntry, PairRecord, FluxEdge, PlaquetteImagResult, SiteFlowResult
export FluxDiagnosticResult, FluxRunResult

struct FluxEndpoint
    value::Union{Int, String}
end

struct FluxRow
    nn::Int
    from::FluxEndpoint
    to::FluxEndpoint
    R::RKey
    value::ComplexF64
end

struct FluxTerm
    rows::Vector{FluxRow}
end

struct RunFiles
    hr_path::String
    win_path::Union{Nothing, String}
    poscar_path::Union{Nothing, String}
end

struct GeometryConfig
    search_bounds::RKey
    distance_tol::Float64
end

struct ArrowStyle
    selector::String
    size::Float64
    color::String
end

struct PlotConfig
    interactive::Bool
    cell_bounds::Union{Nothing, RKey}
    arrow_styles::Vector{ArrowStyle}
end

struct SpinConfig
    layout::Symbol
end

struct BasisConfig
    orbitals_per_atom::Dict{String, Int}
    orbitals_per_species_group::Vector{Int}
end

struct DiagnosticVertex
    endpoint::FluxEndpoint
    cell::RKey
end

struct PlaquetteDiagnostic
    name::String
    vertices::Vector{DiagnosticVertex}
end

struct DiagnosticConfig
    enabled::Bool
    continuity::Bool
    continuity_tol::Float64
    plaquettes::Vector{PlaquetteDiagnostic}
end

struct FluxConfig
    files::RunFiles
    geometry::GeometryConfig
    plot::PlotConfig
    spin::SpinConfig
    basis::BasisConfig
    diagnostic::DiagnosticConfig
    terms::Vector{FluxTerm}
end

struct FluxBasisEntry
    index::Int
    site_label::String
    species::String
    orbital::String
    spin::Symbol
    center_frac::NTuple{3, Float64}
end

struct PairRecord
    i::Int
    j::Int
    R::RKey
    distance::Float64
    shell::Int
end

struct FluxEdge
    nn::Int
    from_index::Int
    to_index::Int
    R::RKey
    value::ComplexF64
    from_label::String
    to_label::String
    start_frac::NTuple{3, Float64}
    finish_frac::NTuple{3, Float64}
    start_cart::NTuple{3, Float64}
    finish_cart::NTuple{3, Float64}
end

struct PlaquetteImagResult
    name::String
    imag_sum::Float64
    edge_imags::Vector{Float64}
end

struct SiteFlowResult
    index::Int
    label::String
    flow_in::Float64
    flow_out::Float64
    residual::Float64
    passed::Bool
end

struct FluxDiagnosticResult
    plaquettes::Vector{PlaquetteImagResult}
    site_flows::Vector{SiteFlowResult}
    continuity_passed::Bool
end

struct FluxRunResult
    config::FluxConfig
    output_hr::String
    html_path::Union{Nothing, String}
    diagnostic_path::Union{Nothing, String}
    edges::Vector{FluxEdge}
    diagnostic::Union{Nothing, FluxDiagnosticResult}
    roundtrip_validated::Bool
end

end
