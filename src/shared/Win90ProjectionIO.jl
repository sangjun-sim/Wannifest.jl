module Win90ProjectionIO

using LinearAlgebra: norm

using ..CrystalCells: wrap_fractional
using ..Win90OrbitalTokens: expand_projection_orbitals

export WinProjectionSeed, WinProjectionSource, read_win_projection_source

const BOHR_TO_ANGSTROM = 0.529177210903

struct WinProjectionSeed
    species::String
    site::Int
    site_label::String
    orbital::String
    center_frac::NTuple{3, Float64}
end

struct WinProjectionSource
    num_wann::Int
    spinors::Bool
    species_atoms::Dict{String, Vector{NTuple{3, Float64}}}
    seeds::Vector{WinProjectionSeed}
    metadata::Dict{String, Any}
end

function _strip_wannier_comment(line::AbstractString)
    bang = findfirst('!', line)
    clean = isnothing(bang) ? String(line) : String(line[begin:prevind(line, bang)])
    hash = findfirst('#', clean)
    clean = isnothing(hash) ? clean : String(clean[begin:prevind(clean, hash)])
    return strip(clean)
end

function _unit_scale(unit::AbstractString, path::AbstractString, line_number::Integer)
    lower = lowercase(strip(unit))
    lower in ("ang", "angstrom", "angstroms") && return 1.0
    lower in ("bohr", "bohrs", "au", "a.u.") && return BOHR_TO_ANGSTROM
    error("Unsupported wannier90 unit '$unit' on line $line_number in $(abspath(path))")
end

function _parse_vec3(tokens, context::AbstractString)
    length(tokens) >= 3 || error("Malformed $context: expected three numeric values")
    return (parse(Float64, tokens[1]), parse(Float64, tokens[2]), parse(Float64, tokens[3]))
end

function _matrix_from_wannier_rows(rows::Vector{NTuple{3, Float64}}, scale::Float64)
    length(rows) == 3 || error("unit_cell_cart must contain exactly three lattice vectors")
    lattice = Matrix{Float64}(undef, 3, 3)
    for j in 1:3
        lattice[:, j] .= scale .* collect(rows[j])
    end
    return lattice
end

function _wrap_frac(v::AbstractVector{<:Real})
    wrapped = wrap_fractional(Float64[v[1], v[2], v[3]])
    return (wrapped[1], wrapped[2], wrapped[3])
end
_tuple_vec(v::NTuple{3, Float64}) = [v[1], v[2], v[3]]

function _frac_distance(a::NTuple{3, Float64}, b::NTuple{3, Float64})
    d = [a[1] - b[1], a[2] - b[2], a[3] - b[3]]
    d .-= round.(d)
    return norm(d)
end

function _parse_bool(raw::AbstractString)
    lower = lowercase(strip(raw))
    lower in ("true", ".true.", "t") && return true
    lower in ("false", ".false.", "f") && return false
    error("spinors must be true/false, got '$raw'")
end

function _assignment_value(clean::AbstractString, key::AbstractString)
    lower = lowercase(clean)
    startswith(lower, lowercase(key)) || return nothing
    occursin("=", clean) || return nothing
    parts = split(clean, '='; limit=2)
    return strip(parts[2])
end

function _site_records(
    lattice::Matrix{Float64},
    atom_labels::Vector{String},
    atom_positions::Vector{NTuple{3, Float64}},
    atom_mode::Symbol,
    atom_scale::Float64,
)
    counters = Dict{String, Int}()
    records = NamedTuple[]
    species_atoms = Dict{String, Vector{NTuple{3, Float64}}}()
    for (label, position) in zip(atom_labels, atom_positions)
        counters[label] = get(counters, label, 0) + 1
        site = counters[label]
        frac = if atom_mode == :cartesian
            _wrap_frac(lattice \ (atom_scale .* _tuple_vec(position)))
        else
            _wrap_frac(_tuple_vec(position))
        end
        push!(get!(species_atoms, label, NTuple{3, Float64}[]), frac)
        push!(records, (
            species=label,
            site=site,
            site_label=string(label, site),
            center_frac=frac,
        ))
    end
    return records, species_atoms
end

function _match_site(records, frac::NTuple{3, Float64}; tol::Float64=1e-3)
    isempty(records) && return nothing
    distances = [_frac_distance(record.center_frac, frac) for record in records]
    idx = argmin(distances)
    return distances[idx] <= tol ? records[idx] : nothing
end

function _parse_projection_lhs(
    lhs::AbstractString,
    lattice::Matrix{Float64},
    site_records,
    unknown_counter::Base.RefValue{Int},
)
    clean = strip(String(lhs))
    lower = lowercase(clean)
    if startswith(lower, "f=") || startswith(lower, "c=")
        mode = startswith(lower, "f=") ? :fractional : :cartesian
        raw_coords = strip(split(clean, '='; limit=2)[2])
        coords = _parse_vec3(split(raw_coords, ','), "projection coordinate")
        frac = mode == :fractional ? _wrap_frac(collect(coords)) : _wrap_frac(lattice \ collect(coords))
        matched = _match_site(site_records, frac)
        if isnothing(matched)
            unknown_counter[] += 1
            return (
                species="unknown",
                site=unknown_counter[],
                site_label=string("unknown", unknown_counter[]),
                center_frac=frac,
            )
        end
        return matched
    end

    matches = [record for record in site_records if record.species == clean]
    isempty(matches) && error("projection species '$clean' is not present in atoms block")
    return matches
end

function _first_projection_colon(line::AbstractString)
    colon = findfirst(==(':'), line)
    isnothing(colon) && error("Malformed projection line without ':' separator: $line")
    lhs = strip(line[begin:prevind(line, colon)])
    rhs = strip(line[nextind(line, colon):end])
    second_colon = findfirst(==(':'), rhs)
    orbitals = isnothing(second_colon) ? rhs : strip(rhs[begin:prevind(rhs, second_colon)])
    return lhs, orbitals
end

function _push_orbital_seeds!(
    out::Vector{WinProjectionSeed},
    site,
    orbitals::Vector{String},
)
    for orbital in orbitals
        push!(out, WinProjectionSeed(site.species, site.site, site.site_label, orbital, site.center_frac))
    end
    return out
end

function read_win_projection_source(path::AbstractString)
    num_wann = nothing
    spinors = false
    lattice_rows = NTuple{3, Float64}[]
    lattice_scale = 1.0
    atom_labels = String[]
    atom_positions = NTuple{3, Float64}[]
    atom_mode = nothing
    atom_scale = 1.0
    projection_lines = String[]
    metadata = Dict{String, Any}("path" => abspath(path), "raw_projection_lines" => projection_lines)

    active_block = nothing
    for (line_number, raw) in enumerate(readlines(path))
        clean = _strip_wannier_comment(raw)
        isempty(clean) && continue
        lower = lowercase(clean)

        if isnothing(active_block)
            if !isnothing(_assignment_value(clean, "num_wann"))
                num_wann = parse(Int, _assignment_value(clean, "num_wann"))
                continue
            elseif !isnothing(_assignment_value(clean, "spinors"))
                spinors = _parse_bool(_assignment_value(clean, "spinors"))
                metadata["spinors_line"] = line_number
                continue
            end

            tokens = split(lower)
            if length(tokens) >= 2 && tokens[1] == "begin"
                block = tokens[2]
                if block == "unit_cell_cart"
                    length(tokens) >= 3 && (lattice_scale = _unit_scale(tokens[3], path, line_number))
                    active_block = :unit_cell_cart
                elseif block == "atoms_cart"
                    atom_mode = :cartesian
                    length(tokens) >= 3 && (atom_scale = _unit_scale(tokens[3], path, line_number))
                    active_block = :atoms_cart
                elseif block == "atoms_frac"
                    atom_mode = :fractional
                    active_block = :atoms_frac
                elseif block == "projections"
                    active_block = :projections
                end
            end
            continue
        end

        if startswith(lower, "end")
            active_block = nothing
            continue
        end

        if active_block == :unit_cell_cart
            if isempty(lattice_rows) && length(split(clean)) == 1
                lattice_scale = _unit_scale(clean, path, line_number)
                continue
            end
            push!(lattice_rows, _parse_vec3(split(clean), "unit_cell_cart vector"))
        elseif active_block == :atoms_cart || active_block == :atoms_frac
            tokens = split(clean)
            length(tokens) >= 4 || error("Malformed atoms block line $line_number in $(abspath(path))")
            push!(atom_labels, String(tokens[1]))
            push!(atom_positions, _parse_vec3(tokens[2:4], "atom position"))
        elseif active_block == :projections
            occursin("random", lowercase(clean)) && error("random projections are not supported by mode=\"win_groups\"")
            push!(projection_lines, clean)
        end
    end

    isnothing(num_wann) && error("Missing num_wann in $(abspath(path))")
    isempty(projection_lines) && error("Missing or empty projections block in $(abspath(path))")
    isempty(lattice_rows) && error("Missing unit_cell_cart block in $(abspath(path))")
    isnothing(atom_mode) && error("Missing atoms_cart or atoms_frac block in $(abspath(path))")
    spinors && isodd(num_wann) && error("spinors=.true. requires even num_wann, got $num_wann")

    lattice = _matrix_from_wannier_rows(lattice_rows, lattice_scale)
    site_records, species_atoms = _site_records(lattice, atom_labels, atom_positions, atom_mode, atom_scale)
    unknown_counter = Ref(0)
    seeds = WinProjectionSeed[]

    for line in projection_lines
        lhs, orbital_text = _first_projection_colon(line)
        projection_orbitals = expand_projection_orbitals(orbital_text)
        sites = _parse_projection_lhs(lhs, lattice, site_records, unknown_counter)
        if sites isa AbstractVector
            for site in sites
                _push_orbital_seeds!(seeds, site, projection_orbitals)
            end
        else
            _push_orbital_seeds!(seeds, sites, projection_orbitals)
        end
    end

    return WinProjectionSource(Int(num_wann), spinors, species_atoms, seeds, metadata)
end

end
