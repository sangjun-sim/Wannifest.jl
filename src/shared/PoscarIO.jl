module PoscarIO

using LinearAlgebra
using Printf

using ..CrystalCells: CrystalCell, natoms, species_counts, species_labels, wrap_fractional

export StructureCell
export natoms
export species_counts
export species_labels
export wrap_fractional
export read_poscar
export write_poscar

const StructureCell = CrystalCell

function scaled_lattice(raw_lattice::Matrix{Float64}, scale_tokens::Vector{SubString{String}})
    if length(scale_tokens) == 1
        scale = parse(Float64, scale_tokens[1])
        if scale > 0
            return scale .* raw_lattice
        elseif scale < 0
            volume = abs(det(raw_lattice))
            volume > 0 || error("Lattice vectors are singular in the structure file")
            factor = cbrt(abs(scale) / volume)
            return factor .* raw_lattice
        end
        error("The POSCAR scaling factor cannot be zero")
    elseif length(scale_tokens) == 3
        scales = parse.(Float64, scale_tokens)
        all(scales .> 0) || error("Three-component POSCAR scaling factors must be positive")
        return raw_lattice * Diagonal(scales)
    end

    error("Unsupported POSCAR scaling specification on line 2")
end

function _is_integer_tokens(tokens::Vector{SubString{String}})
    isempty(tokens) && return false
    try
        parse.(Int, tokens)
        return true
    catch
        return false
    end
end

function _ordered_atom_indices(cell::StructureCell)
    grouped = Int[]
    for species_id in eachindex(cell.species_names)
        for atom_index in eachindex(cell.species_ids)
            if cell.species_ids[atom_index] == species_id
                push!(grouped, atom_index)
            end
        end
    end
    return grouped
end

function _format_triplet(values::AbstractVector{<:Real}, digits::Int)
    return @sprintf("% .*f  % .*f  % .*f", digits, float(values[1]), digits, float(values[2]), digits, float(values[3]))
end

function read_poscar(path::AbstractString)::StructureCell
    lines = readlines(path)
    length(lines) >= 8 || error("Structure file is too short: $(abspath(path))")

    comment = rstrip(lines[1])
    scale_tokens = split(strip(lines[2]))
    raw_lattice = Matrix{Float64}(undef, 3, 3)

    for j in 1:3
        tokens = split(strip(lines[2 + j]))
        length(tokens) >= 3 || error("Malformed lattice vector on line $(2 + j) in $(abspath(path))")
        raw_lattice[:, j] = parse.(Float64, tokens[1:3])
    end

    lattice = scaled_lattice(raw_lattice, scale_tokens)

    index = 6
    first_tokens = split(strip(lines[index]))
    species_names, counts = if _is_integer_tokens(first_tokens)
        (["Type$(i)" for i in eachindex(first_tokens)], parse.(Int, first_tokens))
    else
        parsed_species_names = String.(first_tokens)
        index += 1
        index <= length(lines) || error("Missing atom counts line in $(abspath(path))")
        count_tokens = split(strip(lines[index]))
        _is_integer_tokens(count_tokens) || error("Expected atom counts on line $(index) in $(abspath(path))")
        (parsed_species_names, parse.(Int, count_tokens))
    end

    index += 1
    index <= length(lines) || error("Missing coordinate mode line in $(abspath(path))")
    mode_line = strip(lines[index])
    if !isempty(mode_line) && startswith(lowercase(mode_line), "s")
        error("Selective dynamics is not supported in this project: $(abspath(path))")
    end

    lower_mode = lowercase(mode_line)
    direct_mode = startswith(lower_mode, "d")
    cartesian_mode = startswith(lower_mode, "c") || startswith(lower_mode, "k")
    (direct_mode || cartesian_mode) || error("Unsupported coordinate mode on line $(index) in $(abspath(path))")

    index += 1
    total_atoms = sum(counts)
    index + total_atoms - 1 <= length(lines) || error("Expected $total_atoms atomic positions in $(abspath(path))")

    coordinates = Matrix{Float64}(undef, 3, total_atoms)
    for atom_index in 1:total_atoms
        tokens = split(strip(lines[index + atom_index - 1]))
        length(tokens) >= 3 || error("Malformed coordinate line $(index + atom_index - 1) in $(abspath(path))")
        coordinates[:, atom_index] = parse.(Float64, tokens[1:3])
    end

    frac_positions = if direct_mode
        coordinates
    else
        lattice \ coordinates
    end

    species_ids = Int[]
    for (species_id, count) in enumerate(counts)
        append!(species_ids, fill(species_id, count))
    end

    return StructureCell(
        comment,
        lattice,
        wrap_fractional(frac_positions),
        species_names,
        species_ids,
        abspath(path),
    )
end

function write_poscar(path::AbstractString, cell::StructureCell; digits::Int=12)
    order = _ordered_atom_indices(cell)
    counts = species_counts(cell)
    mkpath(dirname(abspath(path)))

    open(path, "w") do io
        println(io, cell.comment)
        println(io, "1.0")
        for j in 1:3
            println(io, _format_triplet(cell.lattice[:, j], digits))
        end
        println(io, join(cell.species_names, " "))
        println(io, join(string.(counts), " "))
        println(io, "Direct")
        for atom_index in order
            x = wrap_fractional(cell.frac_positions[:, atom_index])
            println(io, _format_triplet(x, digits))
        end
    end

    return path
end

end
