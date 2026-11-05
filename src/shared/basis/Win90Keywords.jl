module Win90Keywords

export strip_wannier_comment, win_keyword_assignment, read_win_keywords, read_win_block
export parse_win_bool, parse_win_int, parse_win_real, parse_band_ranges

strip_wannier_comment(line::AbstractString) = strip(first(split(line, r"[!#]"; limit=2)))

const _ASSIGNMENT = r"^([A-Za-z_][A-Za-z0-9_]*)(?:\s*[=:]\s*|\s+)(.*)$"

"""
    win_keyword_assignment(clean, keys) -> (key, value) | nothing

Match one comment-stripped line against `key = value | key : value | key value`.
Keys are case insensitive, as in Wannier90; the returned key is lowercase.
"""
function win_keyword_assignment(clean::AbstractString, keys)
    matched = match(_ASSIGNMENT, clean)
    isnothing(matched) && return nothing
    key = lowercase(matched.captures[1])
    key in keys || return nothing
    return key, String(strip(matched.captures[2]))
end

function _block_marker(clean::AbstractString, marker::AbstractString)
    tokens = split(lowercase(clean))
    length(tokens) >= 2 && tokens[1] == marker || return nothing
    return String(tokens[2])
end

"""
    read_win_keywords(path, keys) -> Dict{String, String}

Collect the requested keywords found outside `begin ... end` blocks.
A keyword set more than once is an error, as in Wannier90.
"""
function read_win_keywords(path::AbstractString, keys::Set{String})::Dict{String, String}
    isfile(path) || error("Wannier90 input file not found: $(abspath(path))")
    found = Dict{String, String}()
    active_block = nothing
    for raw in eachline(path)
        clean = strip_wannier_comment(raw)
        isempty(clean) && continue
        if isnothing(active_block)
            active_block = _block_marker(clean, "begin")
            isnothing(active_block) || continue
            assignment = win_keyword_assignment(clean, keys)
            isnothing(assignment) && continue
            key, value = assignment
            haskey(found, key) && error("Duplicate $key setting in $(abspath(path))")
            found[key] = value
        elseif _block_marker(clean, "end") == active_block
            active_block = nothing
        end
    end
    isnothing(active_block) || error("Unterminated begin $active_block block in $(abspath(path))")
    return found
end

"""
    read_win_block(path, name) -> Vector{String} | nothing

Return the comment-stripped, non-empty lines of `begin name ... end name`.
"""
function read_win_block(path::AbstractString, name::AbstractString)
    isfile(path) || error("Wannier90 input file not found: $(abspath(path))")
    target = lowercase(name)
    lines = nothing
    active_block = nothing
    for raw in eachline(path)
        clean = strip_wannier_comment(raw)
        isempty(clean) && continue
        if isnothing(active_block)
            active_block = _block_marker(clean, "begin")
            if active_block == target
                isnothing(lines) || error("Duplicate begin $target block in $(abspath(path))")
                lines = String[]
            end
        elseif _block_marker(clean, "end") == active_block
            active_block = nothing
        elseif active_block == target
            push!(lines, String(clean))
        end
    end
    isnothing(active_block) || error("Unterminated begin $active_block block in $(abspath(path))")
    return lines
end

function parse_win_bool(raw::AbstractString)::Bool
    lower = lowercase(strip(raw))
    lower in ("true", ".true.", "t") && return true
    lower in ("false", ".false.", "f") && return false
    error("Invalid Wannier90 logical value '$raw' (expected .true./true/t or .false./false/f)")
end

function parse_win_int(raw::AbstractString)::Int
    value = tryparse(Int, strip(raw))
    isnothing(value) && error("Invalid Wannier90 integer value '$raw'")
    return value
end

"""Parse a real value, accepting Fortran exponents such as `1.0d0` and `-2.5D-1`."""
function parse_win_real(raw::AbstractString)::Float64
    value = tryparse(Float64, replace(strip(raw), r"[dD]" => "e"))
    !isnothing(value) && isfinite(value) || error("Invalid Wannier90 real value '$raw'")
    return value
end

"""
    parse_band_ranges(raw) -> Vector{Int}

Expand `"1-4, 9"`, `"1-4 9"` or `"1:4,9"` into `[1, 2, 3, 4, 9]`.
Entries are separated by commas, semicolons or blanks. Duplicate bands,
reversed ranges and non-positive indices are errors.
"""
function parse_band_ranges(raw::AbstractString)::Vector{Int}
    text = replace(strip(raw), r"\s*[-:]\s*" => "-")
    bands = Int[]
    for token in split(text, r"[\s,;]+"; keepempty=false)
        matched = match(r"^(\d+)(?:-(\d+))?$", token)
        isnothing(matched) && error("Invalid band range entry '$token' in '$raw'")
        first_band = parse(Int, matched.captures[1])
        last_band = isnothing(matched.captures[2]) ? first_band : parse(Int, matched.captures[2])
        first_band >= 1 || error("Band indices must be positive in '$raw'")
        last_band >= first_band || error("Reversed band range '$token' in '$raw'")
        append!(bands, first_band:last_band)
    end
    isempty(bands) && error("Empty band range '$raw'")
    allunique(bands) || error("Duplicate band index in '$raw'")
    return bands
end

end
