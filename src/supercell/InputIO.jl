module InputIO

using LinearAlgebra

export RunConfig
export DEFAULT_INPUT_PATH, DEFAULT_OUTPUT_PATH
export default_config, parse_matrix_args

struct RunConfig
    input_path::String
    output_path::String
    mismatch_output::Union{Nothing, String}
    basis::Symbol
    use_symmetry::Bool
    hall_number::Union{Nothing, Int}
    symprec::Float64
    angle_tolerance::Float64
    digits::Int
    validate::Bool
    matrix_rows::Matrix{Int}
end

const DEFAULT_INPUT_PATH = "POSCAR"
const DEFAULT_OUTPUT_PATH = "POSCAR.supercell"
const DEFAULT_BASIS = :input
const DEFAULT_USE_SYMMETRY = true
const DEFAULT_HALL_NUMBER = nothing
const DEFAULT_SYMPREC = 1.0e-5
const DEFAULT_ANGLE_TOLERANCE = -1.0
const DEFAULT_DIGITS = 12
const DEFAULT_VALIDATE = true

function _reject_extra_args(args::Vector{String})
    isempty(args) && error("Missing supercell matrix. Use either `2 2 1` or `2 0 0, 0 2 0, 0 0 1`.")
    for arg in args
        startswith(arg, "--") && error("supercell accepts only matrix arguments; got \"$arg\"")
        occursin(";", arg) && error("Use commas, not semicolons, between full-matrix rows.")
    end
    return nothing
end

function _parse_int_token(token::AbstractString, context::AbstractString)::Int
    try
        return parse(Int, token)
    catch
        error("$context must contain only integers; got \"$token\"")
    end
end

function _assert_full_rank(matrix::Matrix{Int})::Matrix{Int}
    rank(Matrix{Float64}(matrix)) == 3 || error("supercell matrix must be full-rank")
    return matrix
end

function _diagonal_matrix(values::Vector{Int})::Matrix{Int}
    matrix = zeros(Int, 3, 3)
    for i in 1:3
        matrix[i, i] = values[i]
    end
    return matrix
end

function _parse_diag_matrix(args::Vector{String})::Matrix{Int}
    length(args) == 3 || error("Diagonal supercell form must be exactly three integer entries, e.g. `julia main.jl supercell 2 2 1`.")
    any(arg -> occursin(",", arg), args) && error("Diagonal supercell form must not contain commas.")
    values = [_parse_int_token(token, "Diagonal supercell form") for token in args]
    return _assert_full_rank(_diagonal_matrix(values))
end

function _parse_full_matrix(args::Vector{String})::Matrix{Int}
    length(args) == 9 || error("Full supercell matrix form must be exactly nine entries with commas after entries 3 and 6, e.g. `julia main.jl supercell 2 0 0, 0 2 0, 0 0 1`.")
    endswith(args[3], ",") || error("Full supercell matrix requires a comma after row 1.")
    endswith(args[6], ",") || error("Full supercell matrix requires a comma after row 2.")
    for i in (1, 2, 4, 5, 7, 8, 9)
        occursin(",", args[i]) && error("Full supercell matrix commas are allowed only after entries 3 and 6.")
    end

    matrix = Matrix{Int}(undef, 3, 3)
    tokens = copy(args)
    tokens[3] = chop(tokens[3]; tail=1)
    tokens[6] = chop(tokens[6]; tail=1)
    for i in 1:3, j in 1:3
        matrix[i, j] = _parse_int_token(tokens[3 * (i - 1) + j], "Full supercell matrix")
    end

    return _assert_full_rank(matrix)
end

function parse_matrix_args(args::Vector{String})::Matrix{Int}
    _reject_extra_args(args)
    return any(arg -> occursin(",", arg), args) ? _parse_full_matrix(args) : _parse_diag_matrix(args)
end

function default_config(matrix_rows::Matrix{Int})::RunConfig
    return RunConfig(
        DEFAULT_INPUT_PATH,
        DEFAULT_OUTPUT_PATH,
        nothing,
        DEFAULT_BASIS,
        DEFAULT_USE_SYMMETRY,
        DEFAULT_HALL_NUMBER,
        DEFAULT_SYMPREC,
        DEFAULT_ANGLE_TOLERANCE,
        DEFAULT_DIGITS,
        DEFAULT_VALIDATE,
        _assert_full_rank(matrix_rows),
    )
end

end
