module Validate

using LinearAlgebra
using Printf
using Spglib
using ..PoscarIO: StructureCell, natoms
using ..Symmetry: build_dataset, build_spglib_cell, summarize_dataset

export validate_symmetry_mapping
export write_atom_symm_mismatch

function _apply_symmetry_op(rotation, translation, x::AbstractVector{<:Real})
    y = Matrix{Float64}(rotation) * Float64.(x) + Vector{Float64}(translation)
    return y .- floor.(y)
end

function _best_species_match(cell::StructureCell, x::AbstractVector{<:Real}, species_id::Int)
    best_distance = Inf
    best_index = 0
    best_position = zeros(3)
    for atom_index in 1:natoms(cell)
        cell.species_ids[atom_index] == species_id || continue
        d = cell.frac_positions[:, atom_index] - x
        d .-= round.(d)
        distance = norm(cell.lattice * d)
        if distance < best_distance
            best_distance = distance
            best_index = atom_index
            best_position = copy(cell.frac_positions[:, atom_index])
        end
    end
    return best_distance, best_index, best_position
end

function _det3(r::AbstractMatrix{<:Integer})
    return r[1, 1] * (r[2, 2] * r[3, 3] - r[2, 3] * r[3, 2]) -
           r[1, 2] * (r[2, 1] * r[3, 3] - r[2, 3] * r[3, 1]) +
           r[1, 3] * (r[2, 1] * r[3, 2] - r[2, 2] * r[3, 1])
end

function _rotation_order(rotation::AbstractMatrix{<:Integer}; max_order::Int=12)
    identity_matrix = Matrix{Int}(I, size(rotation, 1), size(rotation, 2))
    current = copy(identity_matrix)
    for order in 1:max_order
        current = current * rotation
        current == identity_matrix && return order
    end
    return nothing
end

function _schonflies_label(rotation)
    r = Matrix{Int}(rotation)
    det_r = _det3(r)
    order = _rotation_order(r)
    trace_r = sum(r[i, i] for i in 1:3)

    if det_r == 1
        isnothing(order) && return "C?"
        order == 1 && return "E"
        return "C$(order)"
    end

    det_r == -1 || return "?"
    r == -Matrix{Int}(I, 3, 3) && return "i"
    order == 2 && trace_r == 1 && return "σ"
    trace_r == -2 && return "S3"
    trace_r == -1 && return "S4"
    trace_r == 0 && return "S6"
    isnothing(order) && return "S?"
    return "S$(order)"
end

function _operation_label(rotation, translation)
    r = vec(Int.(rotation))
    t = Vector{Float64}(translation)
    schonflies = _schonflies_label(rotation)
    rotation_text = @sprintf(
        "[%d,%d,%d;%d,%d,%d;%d,%d,%d]",
        r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8], r[9]
    )
    translation_text = @sprintf("[%.6f,%.6f,%.6f]", t[1], t[2], t[3])
    return string("{", schonflies, " | ", translation_text, "}") # R=", rotation_text, ", t=", translation_text)
end

function write_atom_symm_mismatch(path::AbstractString, rows)
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        @printf(io, "%5s %-4s %9s %9s %9s %-45s %9s %9s %9s %12s\n", "Index", "Atom", "R_x", "R_y", "R_z", "Symmetry operator", "R'_x", "R'_y", "R'_z", "Mismatch")
        for row in rows
            @printf(
                io,
                "%5d %-4s %9.6f %9.6f %9.6f %-45s %9.6f %9.6f %9.6f %12.6e\n",
                row.index,
                row.atom_label,
                row.x0[1], row.x0[2], row.x0[3],
                row.op_label,
                row.x1[1], row.x1[2], row.x1[3],
                row.mismatch,
            )
        end
    end
    return path
end

function validate_symmetry_mapping(
    cell::StructureCell;
    symprec::Float64=1e-5,
    mismatch_path::AbstractString=joinpath(@__DIR__, "atom_symm_mismatch.dat"),
)
    dataset = build_dataset(cell; symprec=symprec)
    summary = summarize_dataset(dataset)
    rotations, translations = Spglib.get_symmetry(build_spglib_cell(cell), symprec)

    rows = NamedTuple[]
    mismatches = Float64[]

    for (op_index, (rotation, translation)) in enumerate(zip(rotations, translations))
        op_label = "op$(op_index): " * _operation_label(rotation, translation)
        for atom_index in 1:natoms(cell)
            x0 = cell.frac_positions[:, atom_index]
            transformed = _apply_symmetry_op(rotation, translation, x0)
            mismatch, matched_index, matched_position = _best_species_match(cell, transformed, cell.species_ids[atom_index])
            atom_label = "$(cell.species_names[cell.species_ids[atom_index]])-$(atom_index)"
            push!(
                rows,
                (
                    index = atom_index,
                    atom_label = atom_label,
                    x0 = copy(x0),
                    op_label = op_label,
                    x1 = copy(transformed),
                    mismatch = mismatch,
                    matched_index = matched_index,
                    matched_position = matched_position,
                ),
            )
            push!(mismatches, mismatch)
        end
    end

    write_atom_symm_mismatch(mismatch_path, rows)

    return (
        summary = summary,
        n_operations = length(rotations),
        n_checks = length(mismatches),
        max_mismatch = maximum(mismatches),
        mean_mismatch = sum(mismatches) / length(mismatches),
        tolerance = symprec,
        mismatch_path = mismatch_path,
        rows = rows,
    )
end

end
