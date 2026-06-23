module ProjectionInputGroups

using ..InputParsing: optional_string, reject_unknown_keys, required_string, required_string_vector
using ..ProjectionModel: ProjectionGroupConfig

export parse_projection_groups, projection_color_groups

const ALLOWED_INDEX_GROUP_KEYS = Set(("label", "indices", "color"))
const ALLOWED_WIN_GROUP_KEYS = Set(("label", "species", "sites", "site_labels", "orbitals", "orbital_shells", "spin", "color"))
const DEFAULT_PROJECTION_COLORS = [
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
]

function _string_value(tbl, key::AbstractString, default::AbstractString; context::AbstractString)
    return optional_string(tbl, key; default=default, context=context)
end

function _int_vector_from_raw(raw, context::AbstractString)
    raw isa AbstractVector || error("$context must be an array of integers")
    values = Int[]
    for (i, item) in enumerate(raw)
        item isa Integer || error("$context must contain only integers (bad entry at index $i)")
        push!(values, Int(item))
    end
    isempty(values) && error("$context cannot be empty")
    return values
end

function _required_int_vector(tbl, key::AbstractString; context::AbstractString)
    haskey(tbl, key) || error("$context.$key is required")
    return _int_vector_from_raw(tbl[key], "$context.$key")
end

function _optional_int_vector(tbl, key::AbstractString; context::AbstractString)
    haskey(tbl, key) || return Int[]
    return _required_int_vector(tbl, key; context=context)
end

function _optional_string_vector(tbl, key::AbstractString; context::AbstractString)
    haskey(tbl, key) || return String[]
    return required_string_vector(tbl, key; context=context)
end

function _string_vector_from_raw(raw, context::AbstractString)
    raw isa AbstractVector || error("$context must be an array of strings")
    values = String[]
    for (i, item) in enumerate(raw)
        item isa AbstractString || error("$context must contain only strings (bad entry at index $i)")
        text = strip(String(item))
        isempty(text) && error("$context cannot contain empty strings")
        push!(values, text)
    end
    isempty(values) && error("$context cannot be empty")
    return values
end

function _projection_group_color(raw, index::Integer, context::AbstractString)
    default = DEFAULT_PROJECTION_COLORS[mod1(index, length(DEFAULT_PROJECTION_COLORS))]
    isnothing(raw) && return default
    raw isa AbstractString || error("$context color must be a string")
    text = strip(String(raw))
    return isempty(text) ? default : text
end

function _allowed_group_keys(mode::Symbol)
    mode == :index_groups && return ALLOWED_INDEX_GROUP_KEYS
    mode == :win_groups && return ALLOWED_WIN_GROUP_KEYS
    error("unsupported projection mode: $mode")
end

function _parse_projection_group_table(tbl, mode::Symbol, index::Integer)::ProjectionGroupConfig
    tbl isa AbstractDict || error("projection group $index must be a table")
    context = "projection.groups[$index]"
    reject_unknown_keys(tbl, _allowed_group_keys(mode), context)

    label = required_string(tbl, "label"; context=context)
    color = _string_value(tbl, "color", DEFAULT_PROJECTION_COLORS[mod1(index, length(DEFAULT_PROJECTION_COLORS))]; context=context)
    indices = mode == :index_groups ? _required_int_vector(tbl, "indices"; context=context) : Int[]
    species = mode == :win_groups ? _optional_string_vector(tbl, "species"; context=context) : String[]
    sites = mode == :win_groups ? _optional_int_vector(tbl, "sites"; context=context) : Int[]
    site_labels = mode == :win_groups ? _optional_string_vector(tbl, "site_labels"; context=context) : String[]
    orbitals = mode == :win_groups ? _optional_string_vector(tbl, "orbitals"; context=context) : String[]
    orbital_shells = mode == :win_groups ? _optional_string_vector(tbl, "orbital_shells"; context=context) : String[]
    spin = mode == :win_groups ? _string_value(tbl, "spin", "any"; context=context) : "any"

    return ProjectionGroupConfig(
        label,
        color,
        indices,
        species,
        sites,
        site_labels,
        orbitals,
        orbital_shells,
        spin,
    )
end

function _parse_projection_group_row(row, mode::Symbol, index::Integer)::ProjectionGroupConfig
    context = "projection.groups[$index]"
    row isa AbstractVector || error("$context must be a table or compact array")
    if mode == :index_groups
        length(row) in (2, 3) || error("$context compact index group must be [label, indices, color?]")
        label = row[1] isa AbstractString ? strip(String(row[1])) : error("$context label must be a string")
        isempty(label) && error("$context label cannot be empty")
        indices = _int_vector_from_raw(row[2], "$context indices")
        color = _projection_group_color(length(row) == 3 ? row[3] : nothing, index, context)
        return ProjectionGroupConfig(label, color, indices, String[], Int[], String[], String[], String[], "any")
    elseif mode == :win_groups
        length(row) in (3, 4) || error("$context compact win group must be [label, species, orbitals, color?]")
        label = row[1] isa AbstractString ? strip(String(row[1])) : error("$context label must be a string")
        isempty(label) && error("$context label cannot be empty")
        species = _string_vector_from_raw(row[2], "$context species")
        orbitals = _string_vector_from_raw(row[3], "$context orbitals")
        color = _projection_group_color(length(row) == 4 ? row[4] : nothing, index, context)
        return ProjectionGroupConfig(label, color, Int[], species, Int[], String[], orbitals, String[], "any")
    end
    error("unsupported projection mode: $mode")
end

function _parse_projection_group(raw, mode::Symbol, index::Integer)::ProjectionGroupConfig
    raw isa AbstractDict && return _parse_projection_group_table(raw, mode, index)
    raw isa AbstractVector && return _parse_projection_group_row(raw, mode, index)
    error("projection group $index must be a table or compact array")
end

function _projection_groups_raw(proj_tbl)
    raw = get(proj_tbl, "groups", [])
    raw isa AbstractVector || error("projection.groups must be an array of tables or compact arrays")
    return raw
end

function parse_projection_groups(proj_tbl, mode::Symbol)::Vector{ProjectionGroupConfig}
    groups = ProjectionGroupConfig[
        _parse_projection_group(group_tbl, mode, i)
        for (i, group_tbl) in enumerate(_projection_groups_raw(proj_tbl))
    ]
    isempty(groups) && error("projection.groups requires at least one group when projection.enabled=true")
    labels = [group.label for group in groups]
    length(unique(labels)) == length(labels) || error("projection.groups contains duplicate labels")
    return groups
end

function projection_color_groups(proj_tbl, labels::Vector{String})::Vector{String}
    if !haskey(proj_tbl, "color_group")
        return copy(labels)
    end

    raw = proj_tbl["color_group"]
    raw isa AbstractVector || error("projection.color_group must be an array of strings")
    groups = required_string_vector(proj_tbl, "color_group"; context="projection")

    length(unique(groups)) == length(groups) || error("projection.color_group contains duplicate labels")
    missing = [group for group in groups if !(group in labels)]
    isempty(missing) ||
        error("projection.color_group contains unknown label(s): $(join(missing, ", "))")
    return groups
end

end
