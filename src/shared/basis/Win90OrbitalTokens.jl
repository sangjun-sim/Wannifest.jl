module Win90OrbitalTokens

export expand_orbital_token, expand_orbitals_and_shells, expand_projection_orbitals

const ORBITAL_ALIASES = Dict{String, Vector{String}}(
    "s" => ["s"],
    "p" => ["pz", "px", "py"],
    "d" => ["dz2", "dxz", "dyz", "dx2-y2", "dxy"],
    "f" => ["fz3", "fxz2", "fyz2", "fz(x2-y2)", "fxyz", "fx(x2-3y2)", "fy(3x2-y2)"],
    "t2g" => ["dxy", "dyz", "dxz"],
    "eg" => ["dz2", "dx2-y2"],
    "sp" => ["sp_1", "sp_2"],
    "sp2" => ["sp2_1", "sp2_2", "sp2_3"],
    "sp3" => ["sp3_1", "sp3_2", "sp3_3", "sp3_4"],
    "sp3d" => ["sp3d_1", "sp3d_2", "sp3d_3", "sp3d_4", "sp3d_5"],
    "sp3d2" => ["sp3d2_1", "sp3d2_2", "sp3d2_3", "sp3d2_4", "sp3d2_5", "sp3d2_6"],
)

const EXACT_ORBITALS = Set(vcat(collect(values(ORBITAL_ALIASES))...))
const VALID_ORBITAL_TOKENS = Set(vcat(collect(keys(ORBITAL_ALIASES)), collect(EXACT_ORBITALS)))

function _canonical_token(token::AbstractString)
    return replace(lowercase(strip(String(token))), " " => "")
end

function _angular_momentum_orbital(token::String)
    pairs = Dict{String, String}()
    for part in split(token, ',')
        kv = split(part, '='; limit=2)
        length(kv) == 2 || continue
        pairs[strip(kv[1])] = strip(kv[2])
    end
    haskey(pairs, "l") || return nothing
    l = parse(Int, pairs["l"])
    mr = parse(Int, get(pairs, "mr", "1"))
    table = if l == 0
        ["s"]
    elseif l == 1
        ["pz", "px", "py"]
    elseif l == 2
        ["dz2", "dxz", "dyz", "dx2-y2", "dxy"]
    elseif l == 3
        ORBITAL_ALIASES["f"]
    else
        error("Unsupported angular momentum projection token '$token'")
    end
    1 <= mr <= length(table) || error("Angular momentum projection token '$token' has mr out of range")
    return table[mr]
end

function expand_orbital_token(token::AbstractString; label::AbstractString="")
    clean = _canonical_token(token)
    isempty(clean) && error("empty orbital token")
    if haskey(ORBITAL_ALIASES, clean)
        return copy(ORBITAL_ALIASES[clean])
    elseif clean in EXACT_ORBITALS
        return [clean]
    elseif startswith(clean, "l=")
        return [_angular_momentum_orbital(clean)]
    end
    prefix = isempty(label) ? "" : "group '$label': "
    error(prefix * "unknown orbital token '$token' (valid: $(join(sort!(collect(VALID_ORBITAL_TOKENS)), ", ")))")
end

function _split_orbital_tokens(raw::AbstractString)
    normalized = replace(String(raw), ';' => ',')
    parts = [strip(part) for part in split(normalized, ',') if !isempty(strip(part))]
    tokens = String[]
    i = 1
    while i <= length(parts)
        part = parts[i]
        if startswith(_canonical_token(part), "l=") && i < length(parts) && startswith(_canonical_token(parts[i + 1]), "mr=")
            push!(tokens, string(part, ",", parts[i + 1]))
            i += 2
        else
            push!(tokens, part)
            i += 1
        end
    end
    return tokens
end

function expand_projection_orbitals(raw::AbstractString)
    orbitals = String[]
    for token in _split_orbital_tokens(raw)
        lowercase(strip(token)) == "random" && error("random projections are not supported by mode=\"win_groups\"")
        append!(orbitals, expand_orbital_token(token))
    end
    isempty(orbitals) && error("projection entry has no orbitals")
    return orbitals
end

function expand_orbitals_and_shells(
    orbitals::Vector{String},
    orbital_shells::Vector{String};
    label::AbstractString="",
)
    expanded = String[]
    for token in orbitals
        append!(expanded, expand_orbital_token(token; label=label))
    end
    for token in orbital_shells
        append!(expanded, expand_orbital_token(token; label=label))
    end
    return unique(expanded)
end

end
