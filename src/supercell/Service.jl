module Service

using LinearAlgebra

using ..InputIO: RunConfig
using ..PoscarIO
using ..Symmetry: SymmetrySummary, build_dataset, summarize_dataset
using ..Transform
using ..Validate

export SupercellRunResult, run, load_input_cell

struct SupercellRunResult
    input_path::String
    output_path::String
    basis::Symbol
    use_symmetry::Bool
    matrix::Matrix{Int}
    multiplicity::Int
    expected_atoms::Int
    actual_atoms::Int
    expected_ratio::Int
    actual_ratio::Float64
    input_summary::Union{Nothing, SymmetrySummary}
    output_summary::Union{Nothing, SymmetrySummary}
    validation::Union{Nothing, NamedTuple}
end

function _volume(cell::PoscarIO.StructureCell)
    return abs(det(cell.lattice))
end

function load_input_cell(path::AbstractString)::PoscarIO.StructureCell
    if lowercase(splitext(path)[2]) == ".toml"
        error("TOML structure files are not supported in supercell. Use POSCAR or CONTCAR.")
    end
    return PoscarIO.read_poscar(path)
end

function run(config::RunConfig)::SupercellRunResult
    config.digits >= 1 || error("supercell digits must be positive")

    cell = load_input_cell(config.input_path)
    input_summary = nothing
    basis_cell = cell
    if config.use_symmetry
        input_dataset = build_dataset(
            cell;
            symprec=config.symprec,
            angle_tolerance=config.angle_tolerance,
            hall_number=config.hall_number,
        )
        input_summary = summarize_dataset(input_dataset)
    end

    supercell = Transform.build_supercell(basis_cell, config.matrix_rows)
    PoscarIO.write_poscar(config.output_path, supercell; digits=config.digits)

    output_summary = if config.use_symmetry
        output_dataset = build_dataset(supercell; symprec=config.symprec, angle_tolerance=config.angle_tolerance)
        summarize_dataset(output_dataset)
    else
        nothing
    end

    multiplicity = Transform.supercell_multiplicity(config.matrix_rows)
    expected_atoms = PoscarIO.natoms(basis_cell) * multiplicity
    actual_atoms = PoscarIO.natoms(supercell)
    expected_ratio = multiplicity
    actual_ratio = _volume(supercell) / _volume(basis_cell)

    validation = if config.validate
        mismatch_path = isnothing(config.mismatch_output) ? joinpath(dirname(abspath(config.output_path)), "atom_symm_mismatch.dat") : config.mismatch_output
        Validate.validate_symmetry_mapping(supercell; symprec=config.symprec, mismatch_path=mismatch_path)
    else
        nothing
    end

    return SupercellRunResult(
        config.input_path,
        config.output_path,
        config.basis,
        config.use_symmetry,
        config.matrix_rows,
        multiplicity,
        expected_atoms,
        actual_atoms,
        expected_ratio,
        actual_ratio,
        input_summary,
        output_summary,
        validation,
    )
end

end
