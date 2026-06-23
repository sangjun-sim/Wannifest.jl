module ObservablesInput

using ..InputParsing: reject_unknown_keys, required_bool, required_string_vector
using ..ObservablesModel: OamConfig, OamOrbitalSelection, SamConfig
using ..ObservablesModel: disabled_oam_config, disabled_sam_config

export parse_oam_config, parse_sam_config

const ALLOWED_OAM_KEYS = Set(("enabled", "orbitals", "degeneracy_tol", "plot_components"))
const ALLOWED_SAM_KEYS = Set(("enabled", "degeneracy_tol", "plot_components"))
const ALLOWED_OAM_SHELLS = Set(("s", "p", "d", "t2g", "eg"))
const OAM_COMPONENTS = Dict(
    "lx" => :Lx,
    "ly" => :Ly,
    "lz" => :Lz,
    "l_norm" => :L_norm,
    "l2" => :L2,
)
const SAM_COMPONENTS = Dict(
    "sx" => :Sx,
    "sy" => :Sy,
    "sz" => :Sz,
    "s_norm" => :S_norm,
    "s2" => :S2,
)

function _bool_value(tbl, key::AbstractString, default::Bool; context::AbstractString)
    return haskey(tbl, key) ? required_bool(tbl, key; context=context) : default
end

function _positive_float(tbl, key::AbstractString, default::Real; context::AbstractString)
    value = get(tbl, key, default)
    value isa Real || error("$context.$key must be numeric")
    result = Float64(value)
    result > 0 || error("$context.$key must be positive")
    return result
end

function _selection_from_row(row, index::Integer)::OamOrbitalSelection
    context = "oam.orbitals[$index]"
    row isa AbstractVector || error("$context must be [site_label, orbital_shell]")
    length(row) == 2 || error("$context must have exactly 2 items: [site_label, orbital_shell]")
    row[1] isa AbstractString || error("$context site_label must be a string")
    row[2] isa AbstractString || error("$context orbital_shell must be a string")
    site = strip(String(row[1]))
    shell = strip(String(row[2]))
    isempty(site) && error("$context site_label cannot be empty")
    isempty(shell) && error("$context orbital_shell cannot be empty")
    lowercase(shell) in ALLOWED_OAM_SHELLS ||
        error("$context orbital_shell must be one of s, p, d, t2g, eg")
    return OamOrbitalSelection(site, shell)
end

function _parse_selections(tbl)::Vector{OamOrbitalSelection}
    haskey(tbl, "orbitals") || error("oam.orbitals is required when oam.enabled=true")
    raw = tbl["orbitals"]
    raw isa AbstractVector || error("oam.orbitals must be an array of [site_label, orbital_shell] rows")
    selections = OamOrbitalSelection[_selection_from_row(row, i) for (i, row) in enumerate(raw)]
    isempty(selections) && error("oam.orbitals cannot be empty when oam.enabled=true")
    keys = Set{Tuple{String, String}}()
    for selection in selections
        key = (selection.site, lowercase(selection.orbital_shell))
        key in keys && error("oam.orbitals contains duplicate selection $(selection.site):$(selection.orbital_shell)")
        push!(keys, key)
    end
    return selections
end

function _parse_plot_components(tbl)::Vector{Symbol}
    haskey(tbl, "plot_components") || return [:Lz]
    raw = required_string_vector(tbl, "plot_components"; context="oam")
    isempty(raw) && error("oam.plot_components cannot be empty")
    components = Symbol[]
    seen = Set{Symbol}()
    for component in raw
        key = lowercase(strip(String(component)))
        parsed = get(OAM_COMPONENTS, key, nothing)
        isnothing(parsed) && error("oam.plot_components entries must be one of Lx, Ly, Lz, L_norm, L2")
        parsed in seen && error("oam.plot_components contains duplicate component $parsed")
        push!(seen, parsed)
        push!(components, parsed)
    end
    return components
end

function _parse_sam_plot_components(tbl)::Vector{Symbol}
    haskey(tbl, "plot_components") || return [:Sz]
    raw = required_string_vector(tbl, "plot_components"; context="sam")
    isempty(raw) && error("sam.plot_components cannot be empty")
    components = Symbol[]
    seen = Set{Symbol}()
    for component in raw
        key = lowercase(strip(String(component)))
        parsed = get(SAM_COMPONENTS, key, nothing)
        isnothing(parsed) && error("sam.plot_components entries must be one of Sx, Sy, Sz, S_norm, S2")
        parsed in seen && error("sam.plot_components contains duplicate component $parsed")
        push!(seen, parsed)
        push!(components, parsed)
    end
    return components
end

function parse_oam_config(oam_tbl)::OamConfig
    reject_unknown_keys(oam_tbl, ALLOWED_OAM_KEYS, "oam")
    enabled = _bool_value(oam_tbl, "enabled", false; context="oam")
    enabled || return disabled_oam_config()
    degeneracy_tol = _positive_float(oam_tbl, "degeneracy_tol", 1.0e-4; context="oam")
    return OamConfig(true, _parse_selections(oam_tbl), degeneracy_tol, _parse_plot_components(oam_tbl))
end

function parse_sam_config(sam_tbl)::SamConfig
    reject_unknown_keys(sam_tbl, ALLOWED_SAM_KEYS, "sam")
    enabled = _bool_value(sam_tbl, "enabled", false; context="sam")
    enabled || return disabled_sam_config()
    degeneracy_tol = _positive_float(sam_tbl, "degeneracy_tol", 1.0e-4; context="sam")
    return SamConfig(true, degeneracy_tol, _parse_sam_plot_components(sam_tbl))
end

end
