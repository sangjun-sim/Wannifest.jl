module InputParsing

export namespaced_root, required_table, optional_table
export required_string, optional_string, required_string_vector
export optional_float, optional_int, required_bool
export resolve_path, optional_existing_path
export reject_unknown_keys, reject_unknown_sibling_tables
export parse_vec3_float, parse_matrix3_int, parse_int_pair, parse_float_pair

_prefix(context::AbstractString) = isempty(context) ? "" : string(context, ": ")

function namespaced_root(cfg, namespace::AbstractString)
    root = haskey(cfg, namespace) ? cfg[namespace] : cfg
    root isa AbstractDict || error("[$namespace] must be a table")
    return root
end

function reject_unknown_sibling_tables(root, allowed::Set{String}, namespace::AbstractString)
    unknown = String[]
    for (key, value) in root
        key_text = String(key)
        is_table = value isa AbstractDict ||
            (value isa AbstractVector && !isempty(value) && all(item -> item isa AbstractDict, value))
        is_table && !(key_text in allowed) && push!(unknown, key_text)
    end
    sort!(unknown)
    isempty(unknown) || error("Unsupported $namespace table(s): $(join(unknown, ", "))")
    return nothing
end

function reject_unknown_keys(tbl, allowed::Set{String}, context::AbstractString)
    unknown = sort!(setdiff(collect(keys(tbl)), collect(allowed)))
    isempty(unknown) || error("Unsupported $context option(s): $(join(unknown, ", "))")
    return nothing
end

function required_table(tbl, key::AbstractString; context::AbstractString="")
    value = get(tbl, key, nothing)
    value isa AbstractDict || error(_prefix(context) * "Missing [$key] table")
    return value
end

function optional_table(tbl, key::AbstractString; context::AbstractString="")
    value = get(tbl, key, nothing)
    isnothing(value) && return nothing
    value isa AbstractDict || error(_prefix(context) * "Key \"$key\" must be a table")
    return value
end

function required_string(tbl, key::AbstractString; context::AbstractString="")
    haskey(tbl, key) || error(_prefix(context) * "Missing required key \"$key\"")
    value = tbl[key]
    value isa AbstractString || error(_prefix(context) * "Key \"$key\" must be a string")
    text = strip(String(value))
    isempty(text) && error(_prefix(context) * "Key \"$key\" cannot be empty")
    return text
end

function optional_string(tbl, key::AbstractString; default=nothing, context::AbstractString="")
    haskey(tbl, key) || return default
    value = tbl[key]
    value isa AbstractString || error(_prefix(context) * "Key \"$key\" must be a string")
    text = strip(String(value))
    return isempty(text) ? default : text
end

function required_string_vector(tbl, key::AbstractString; context::AbstractString="")
    haskey(tbl, key) || error(_prefix(context) * "Missing required key \"$key\"")
    raw = tbl[key]
    raw isa AbstractVector || error(_prefix(context) * "Key \"$key\" must be an array of strings")
    values = String[]
    for (i, item) in enumerate(raw)
        item isa AbstractString || error(_prefix(context) * "Key \"$key\" must contain only strings (bad entry at index $i)")
        text = strip(String(item))
        isempty(text) && error(_prefix(context) * "Key \"$key\" cannot contain empty strings")
        push!(values, text)
    end
    isempty(values) && error(_prefix(context) * "Key \"$key\" cannot be an empty array")
    return values
end

function optional_float(tbl, key::AbstractString; context::AbstractString="")
    haskey(tbl, key) || return nothing
    value = tbl[key]
    if value isa AbstractString
        text = strip(String(value))
        return isempty(text) ? nothing : parse(Float64, text)
    elseif value isa Real
        return Float64(value)
    end
    error(_prefix(context) * "Key \"$key\" must be a number or empty string")
end

function optional_int(tbl, key::AbstractString; context::AbstractString="")
    haskey(tbl, key) || return nothing
    value = tbl[key]
    value isa Integer || error(_prefix(context) * "Key \"$key\" must be an integer")
    return Int(value)
end

function required_bool(tbl, key::AbstractString; context::AbstractString="")
    haskey(tbl, key) || error(_prefix(context) * "Missing required key \"$key\"")
    value = tbl[key]
    value isa Bool || error(_prefix(context) * "Key \"$key\" must be a boolean")
    return value
end

function resolve_path(base_dir::AbstractString, value::Union{Nothing, AbstractString}; empty_value=nothing)
    isnothing(value) && return empty_value
    text = strip(String(value))
    isempty(text) && return empty_value
    return isabspath(text) ? text : normpath(joinpath(base_dir, text))
end

function optional_existing_path(base_dir::AbstractString, tbl, key::AbstractString; context::AbstractString)
    raw = optional_string(tbl, key; default=nothing, context=context)
    path = resolve_path(base_dir, raw; empty_value=nothing)
    isnothing(path) && return nothing
    isfile(path) || error("$context.$key points to missing file: $path")
    return path
end

function parse_vec3_float(raw, key::AbstractString)::NTuple{3, Float64}
    raw isa AbstractVector || error("$key must be a 3-element array")
    length(raw) == 3 || error("$key must have length 3")
    return (Float64(raw[1]), Float64(raw[2]), Float64(raw[3]))
end

function parse_matrix3_int(raw, key::AbstractString)::Matrix{Int}
    raw isa AbstractVector || error("$key must be an array")
    if length(raw) == 9
        vals = Int.(raw)
        return Matrix{Int}(permutedims(reshape(vals, 3, 3)))
    elseif length(raw) == 3 && all(x -> x isa AbstractVector && length(x) == 3, raw)
        mat = Matrix{Int}(undef, 3, 3)
        for i in 1:3, j in 1:3
            mat[i, j] = Int(raw[i][j])
        end
        return mat
    end
    error("$key must be either a flat 9-element array or a 3x3 nested array")
end

function parse_int_pair(raw, key::AbstractString)::Tuple{Int, Int}
    raw isa AbstractVector || error("$key must be a 2-element array")
    length(raw) == 2 || error("$key must have length 2")
    return (Int(raw[1]), Int(raw[2]))
end

function parse_float_pair(raw, key::AbstractString)::Tuple{Float64, Float64}
    raw isa AbstractVector || error("$key must be a 2-element array")
    length(raw) == 2 || error("$key must have length 2")
    return (Float64(raw[1]), Float64(raw[2]))
end

end
