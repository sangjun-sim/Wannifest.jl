module Diagnostics

using Printf

using ..BasisSource: entry_label
using ..Model: DiagnosticConfig, FluxBasisEntry, FluxDiagnosticResult, FluxEdge
using ..Model: FluxEndpoint, PlaquetteDiagnostic, PlaquetteImagResult, SiteFlowResult
using ..WannierTypes: RKey

export run_diagnostics, write_diagnostic_tsv

_reverse_key(R::RKey)::RKey = (-R[1], -R[2], -R[3])

function _edge_key(i::Int, j::Int, R::RKey)
    return (i, j, R[1], R[2], R[3])
end

function _edge_index(edges::Vector{FluxEdge})
    index = Dict{Tuple{Int, Int, Int, Int, Int}, FluxEdge}()
    for edge in edges
        key = _edge_key(edge.from_index, edge.to_index, edge.R)
        haskey(index, key) && error("duplicate flux diagnostic edge i=$(edge.from_index), j=$(edge.to_index), R=$(edge.R)")
        index[key] = edge
    end
    return index
end

function _selector_indices(label::AbstractString, basis::Vector{FluxBasisEntry})
    clean = strip(String(label))
    indices = [entry.index for entry in basis if entry_label(entry) == clean]
    isempty(indices) || return indices
    indices = [entry.index for entry in basis if entry.site_label == clean]
    isempty(indices) || return indices
    indices = [entry.index for entry in basis if entry.species == clean]
    isempty(indices) || return indices
    error("No Wannier orbital found for diagnostic selector '$clean'")
end

function _resolve_single_index(endpoint::FluxEndpoint, basis::Vector{FluxBasisEntry})
    if endpoint.value isa Int
        idx = endpoint.value::Int
        1 <= idx <= length(basis) || error("Diagnostic Wannier index $idx is outside 1:$(length(basis))")
        return idx
    end
    label = endpoint.value::String
    indices = _selector_indices(label, basis)
    length(indices) == 1 ||
        error("Diagnostic selector '$label' matched multiple orbitals: $(join(sort(indices), ", "))")
    return only(indices)
end

function _sub_R(b::RKey, a::RKey)::RKey
    return (b[1] - a[1], b[2] - a[2], b[3] - a[3])
end

function _closed_pairs(vertices)
    return ((vertices[i], vertices[i == length(vertices) ? 1 : i + 1]) for i in eachindex(vertices))
end

function _plaquette_imag_sum(edge_index, basis::Vector{FluxBasisEntry}, plaquette::PlaquetteDiagnostic)
    edge_imags = Float64[]
    for (a, b) in _closed_pairs(plaquette.vertices)
        i = _resolve_single_index(a.endpoint, basis)
        j = _resolve_single_index(b.endpoint, basis)
        R = _sub_R(b.cell, a.cell)
        key = _edge_key(i, j, R)
        if haskey(edge_index, key)
            push!(edge_imags, imag(edge_index[key].value))
            continue
        end
        reverse_key = _edge_key(j, i, _reverse_key(R))
        haskey(edge_index, reverse_key) ||
            error("missing flux seed for plaquette '$(plaquette.name)' edge $i -> $j, R=$R")
        push!(edge_imags, -imag(edge_index[reverse_key].value))
    end
    return PlaquetteImagResult(plaquette.name, sum(edge_imags), edge_imags)
end

function _basis_labels(basis::Vector{FluxBasisEntry})
    return Dict(entry.index => entry_label(entry) for entry in basis)
end

function _flow_results(edges::Vector{FluxEdge}, basis::Vector{FluxBasisEntry}, tol::Float64)
    labels = _basis_labels(basis)
    flow_in = Dict(entry.index => 0.0 for entry in basis)
    flow_out = Dict(entry.index => 0.0 for entry in basis)
    for edge in edges
        im = imag(edge.value)
        iszero(im) && continue
        amount = abs(im)
        if im > 0
            flow_out[edge.to_index] += amount
            flow_in[edge.from_index] += amount
        else
            flow_out[edge.from_index] += amount
            flow_in[edge.to_index] += amount
        end
    end
    return [
        begin
            residual = flow_in[entry.index] - flow_out[entry.index]
            SiteFlowResult(entry.index, labels[entry.index], flow_in[entry.index], flow_out[entry.index], residual, abs(residual) <= tol)
        end for entry in sort(basis; by=entry -> entry.index)
    ]
end

function run_diagnostics(edges::Vector{FluxEdge}, basis::Vector{FluxBasisEntry}, config::DiagnosticConfig)
    edge_index = _edge_index(edges)
    plaquettes = [_plaquette_imag_sum(edge_index, basis, plaquette) for plaquette in config.plaquettes]
    site_flows = config.continuity ? _flow_results(edges, basis, config.continuity_tol) : SiteFlowResult[]
    continuity_passed = isempty(site_flows) || all(row -> row.passed, site_flows)
    return FluxDiagnosticResult(plaquettes, site_flows, continuity_passed)
end

_fmt(value::Real) = @sprintf("%.12f", Float64(value))

function write_diagnostic_tsv(path::AbstractString, result::FluxDiagnosticResult)
    mkpath(dirname(path))
    open(path, "w") do io
        if !isempty(result.plaquettes)
            println(io, "# plaquette_imag")
            println(io, "name\timag_sum\tedge_count\tedge_imags")
            for row in result.plaquettes
                println(io, row.name, "\t", _fmt(row.imag_sum), "\t", length(row.edge_imags), "\t", join(_fmt.(row.edge_imags), ","))
            end
        end
        if !isempty(result.site_flows)
            isempty(result.plaquettes) || println(io)
            println(io, "# continuity")
            println(io, "index\tlabel\tflow_in\tflow_out\tresidual\tpassed")
            for row in result.site_flows
                println(io, row.index, "\t", row.label, "\t", _fmt(row.flow_in), "\t", _fmt(row.flow_out), "\t", _fmt(row.residual), "\t", row.passed)
            end
        end
    end
    return path
end

end
