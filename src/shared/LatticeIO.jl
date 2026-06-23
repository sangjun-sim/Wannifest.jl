module LatticeIO

using ..CrystalCells: reciprocal_lattice
import ..PoscarIO

export LatticeData, read_lattice, read_poscar, read_wannier_win

const BOHR_TO_ANGSTROM = 0.529177210903

struct LatticeData
    real_lattice::Matrix{Float64}
    reciprocal_lattice::Matrix{Float64}
    source::String
end

function _lattice_data(lattice::AbstractMatrix{<:Real}, source::AbstractString)
    lat = Matrix{Float64}(lattice)
    return LatticeData(lat, reciprocal_lattice(lat), String(source))
end

function _strip_wannier_comment(line::AbstractString)
    bang = findfirst('!', line)
    clean = if isnothing(bang)
        String(line)
    elseif bang == firstindex(line)
        ""
    else
        String(line[begin:prevind(line, bang)])
    end
    hash = findfirst('#', clean)
    if isnothing(hash)
        return strip(clean)
    elseif hash == firstindex(clean)
        return ""
    end
    return strip(clean[begin:prevind(clean, hash)])
end

function _unit_cell_cart_scale(unit::AbstractString, path::AbstractString, line_number::Integer)
    lower = lowercase(strip(unit))
    if lower in ("ang", "angstrom", "angstroms")
        return 1.0
    elseif lower in ("bohr", "bohrs", "au", "a.u.")
        return BOHR_TO_ANGSTROM
    end
    error("Unsupported unit_cell_cart unit '$unit' on line $line_number in $(abspath(path))")
end

function _parse_wannier_lattice_vector(line::AbstractString, path::AbstractString, line_number::Integer)
    tokens = split(line)
    length(tokens) >= 3 || error("Malformed unit_cell_cart vector on line $line_number in $(abspath(path))")
    return parse.(Float64, tokens[1:3])
end

function _lattice_from_wannier_rows(rows::Vector{Vector{Float64}}, scale::Float64)
    lat = Matrix{Float64}(undef, 3, 3)
    for j in 1:3
        lat[:, j] .= scale .* rows[j]
    end
    return lat
end

function read_wannier_win(path::AbstractString)
    lattice_rows = Vector{Float64}[]
    lattice_scale = 1.0
    active_block = nothing
    found_lattice = false
    found_lattice_end = false

    for (line_number, raw) in enumerate(readlines(path))
        clean = _strip_wannier_comment(raw)
        isempty(clean) && continue
        lower = lowercase(clean)

        if isnothing(active_block)
            tokens = split(lower)
            if length(tokens) >= 2 && tokens[1] == "begin"
                length(tokens) <= 3 || error("Malformed unit_cell_cart begin line $line_number in $(abspath(path))")
                block = tokens[2]
                if block == "unit_cell_cart"
                    if length(tokens) == 3
                        lattice_scale = _unit_cell_cart_scale(tokens[3], path, line_number)
                    end
                    active_block = :unit_cell_cart
                    found_lattice = true
                end
            end
            continue
        end

        if startswith(lower, "end")
            tokens = split(lower)
            length(tokens) >= 2 && tokens[1] == "end" || error("Unexpected end marker on line $line_number in $(abspath(path)): $clean")
            expected = active_block === :unit_cell_cart ? "unit_cell_cart" :
                error("Unexpected active block $active_block")
            tokens[2] == expected || error("Unexpected end marker on line $line_number in $(abspath(path)): $clean")
            active_block === :unit_cell_cart && (found_lattice_end = true)
            active_block = nothing
            continue
        end

        if active_block === :unit_cell_cart
            if isempty(lattice_rows) && length(split(clean)) == 1
                lattice_scale = _unit_cell_cart_scale(clean, path, line_number)
                continue
            end

            length(lattice_rows) < 3 || error("unit_cell_cart has more than three lattice vectors in $(abspath(path))")
            push!(lattice_rows, _parse_wannier_lattice_vector(clean, path, line_number))
            continue
        end
    end

    found_lattice || error("Missing begin unit_cell_cart block in $(abspath(path))")
    found_lattice_end || error("Missing end unit_cell_cart block in $(abspath(path))")
    length(lattice_rows) == 3 || error("unit_cell_cart must contain exactly three lattice vectors in $(abspath(path))")

    lat = _lattice_from_wannier_rows(lattice_rows, lattice_scale)
    return _lattice_data(lat, string(abspath(path), "::unit_cell_cart"))
end

function read_poscar(path::AbstractString)
    cell = PoscarIO.read_poscar(path)
    return _lattice_data(cell.lattice, abspath(path))
end

function read_lattice(path::AbstractString)
    ext = lowercase(splitext(path)[2])
    ext == ".toml" && error("TOML structure files are not supported for band lattice input. Use POSCAR/CONTCAR or wannier90.win.")
    ext == ".win" && return read_wannier_win(path)
    return read_poscar(path)
end

end
