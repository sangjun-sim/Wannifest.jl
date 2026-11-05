module QeInputCheck

export system_tokens, check_no_symmetry, tilted_angles, warn_tilted_magnetization, check_qe_inputs

const TOKEN = r"""![^\r\n]*|'(?:[^']|'')*'|"(?:[^"]|"")*"|[&/=,()%]|[^\s!'"&/=,()%]+|['"]"""
const TRUE_VALUES = (".true.", ".t.", "true", "t")

"""
    system_tokens(text, path) -> Vector{String}

Tokens of the single `&SYSTEM` namelist of a QE input. `!` comments are dropped,
quoted strings stay opaque and names are lowercased; `path` labels error messages.
"""
function system_tokens(text::AbstractString, path::AbstractString)
    tokens = String[]
    for matched in eachmatch(TOKEN, text)
        token = matched.match
        startswith(token, "!") && continue
        if startswith(token, "'") || startswith(token, "\"")
            length(token) > 1 || error("$path: unterminated quoted string in QE input")
            push!(tokens, "<quoted-string>")
        else
            push!(tokens, lowercase(token))
        end
    end
    starts = [i for i in 1:(length(tokens) - 1) if tokens[i] == "&" && tokens[i + 1] == "system"]
    length(starts) == 1 || error("$path: expected exactly one &SYSTEM namelist containing nosym and noinv")
    first = only(starts) + 2
    for i in first:length(tokens)
        if tokens[i] == "/" || (tokens[i] == "&" && i < length(tokens) && tokens[i + 1] == "end")
            return tokens[first:(i - 1)]
        elseif tokens[i] == "&"
            error("$path: &SYSTEM must end before the next namelist")
        end
    end
    error("$path: unterminated &SYSTEM namelist; expected / or &END")
end

function _read_system_tokens(input_file::AbstractString)
    isfile(input_file) || error("$input_file: QE input file not found; nosym and noinv cannot be checked")
    return system_tokens(read(input_file, String), input_file)
end

function _check_symmetry_flags(tokens, input_file)
    for flag in ("nosym", "noinv")
        positions = findall(==(flag), tokens)
        length(positions) == 1 ||
            error("$input_file: &SYSTEM must explicitly set $flag=.true. exactly once (found $(length(positions)))")
        i = only(positions)
        i + 2 <= length(tokens) && tokens[i + 1] == "=" ||
            error("$input_file: malformed $flag assignment; require $flag=.true.")
        tokens[i + 2] in TRUE_VALUES || error("$input_file: $flag must be true; found $(tokens[i + 2])")
        if i + 3 <= length(tokens)
            next_assignment = i + 4 <= length(tokens) &&
                occursin(r"^[a-z_][a-z0-9_]*$", tokens[i + 3]) && tokens[i + 4] in ("=", "(")
            (tokens[i + 3] == "," || next_assignment) ||
                error("$input_file: invalid value after $flag; require a single true value")
        end
    end
    return nothing
end

"""
    check_no_symmetry(input_file)

Require one explicit true assignment each for `nosym` and `noinv` in `&SYSTEM`,
so that QE wrote the full k list without symmetry reduction. This checks the
supplied input text, not its provenance.
"""
check_no_symmetry(input_file::AbstractString) = _check_symmetry_flags(_read_system_tokens(input_file), input_file)

"""Nonzero `angle1(i)` / `angle2(i)` entries of `&SYSTEM`, as `"angle1(2)" => value` pairs."""
function tilted_angles(tokens::Vector{String})
    found = Pair{String, Float64}[]
    for i in 1:(length(tokens) - 5)
        tokens[i] in ("angle1", "angle2") && tokens[i + 1] == "(" && tokens[i + 3] == ")" &&
            tokens[i + 4] == "=" || continue
        value = tryparse(Float64, replace(tokens[i + 5], r"[dD]" => "e"))
        isnothing(value) && continue
        iszero(value) || push!(found, "$(tokens[i])($(tokens[i + 2]))" => value)
    end
    return found
end

"""
    warn_tilted_magnetization(soc_input) -> Int

V pairs the SP up/down states with spinor components 1/2, which are quantised
along z. Warn when the SOC input starts the magnetization away from z; return the
number of nonzero angles.
"""
function warn_tilted_magnetization(soc_input::AbstractString)::Int
    found = tilted_angles(_read_system_tokens(soc_input))
    isempty(found) && return 0
    listing = join(("$name = $value" for (name, value) in found), ", ")
    @warn "SOC input sets $listing; V assumes the SOC magnetization along z, so SP up/down no longer match spinor components 1/2." input=soc_input
    return length(found)
end

"""Check both QE inputs before any calculation file is opened."""
function check_qe_inputs(sp_input::AbstractString, soc_input::AbstractString)
    check_no_symmetry(sp_input)
    check_no_symmetry(soc_input)
    return warn_tilted_magnetization(soc_input)
end

end
