module WannierUmatIO

export UMatrix, read_umat_file

struct UMatrix
    kpoints::Matrix{Float64}        # 3×nk, fractional reciprocal coordinates
    data::Array{ComplexF64, 3}      # nrow×ncol×nk, exactly as stored in the file
    source::String
end

function _next_fields(io::IO, path::AbstractString, label::AbstractString)
    while !eof(io)
        line = strip(readline(io))
        isempty(line) || return split(line)
    end
    error("$path: missing $label (truncated U matrix file)")
end

function _reals(fields, count::Int, path::AbstractString, label::AbstractString)
    length(fields) == count || error("$path: $label must contain $count numbers")
    return map(fields) do field
        value = tryparse(Float64, replace(field, r"[dD]" => "e"))
        !isnothing(value) && isfinite(value) || error("$path: invalid $label value $field")
        value
    end
end

function _dimensions(fields, path::AbstractString)
    length(fields) == 3 || error("$path: expected three integer dimensions nk ncol nrow")
    dims = tryparse.(Int, fields)
    all(x -> !isnothing(x) && x > 0, dims) || error("$path: invalid U matrix dimensions")
    return Int.(dims)
end

function _check_expected(actual::Int, expected, name::AbstractString, path::AbstractString)
    isnothing(expected) || actual == expected ||
        error("$path: $name=$actual does not match the expected $expected")
    return nothing
end

"""
    read_umat_file(path; expected_rows=nothing, expected_cols=nothing) -> UMatrix

Read a Wannier90 `<seed>_u.mat` or `<seed>_u_dis.mat` file. The header holds
`nk ncol nrow`; within each k block the row index varies fastest, one complex
number per line. `data[m, n, k]` and the k order are kept as stored: no
conjugation, transposition, normalization or row realignment is applied.
"""
function read_umat_file(path::AbstractString; expected_rows=nothing, expected_cols=nothing)::UMatrix
    isfile(path) || error("U matrix file not found: $(abspath(path)); run Wannier90 with write_u_matrices = .true.")
    open(path, "r") do io
        eof(io) && error("$path: empty U matrix file")
        readline(io) # Creation comment, not numerical data.
        nk, ncol, nrow = _dimensions(_next_fields(io, path, "nk ncol nrow header"), path)
        nrow >= ncol || error("$path: a U matrix needs at least as many rows ($nrow) as columns ($ncol)")
        _check_expected(nrow, expected_rows, "rows", path)
        _check_expected(ncol, expected_cols, "columns", path)
        # Even minimally formatted numeric rows need four bytes per complex entry.
        big(nk) * (5 + 4 * big(nrow) * ncol) <= filesize(io) - position(io) ||
            error("$path: dimensions exceed the available U matrix data")
        kpoints = Matrix{Float64}(undef, 3, nk)
        data = Array{ComplexF64}(undef, nrow, ncol, nk)
        for k in 1:nk
            kpoints[:, k] = _reals(_next_fields(io, path, "k-point $k"), 3, path, "k-point $k")
            for column in 1:ncol, row in 1:nrow
                label = "U[$row,$column,$k]"
                re, im = _reals(_next_fields(io, path, label), 2, path, label)
                data[row, column, k] = complex(re, im)
            end
        end
        all(line -> isempty(strip(line)), eachline(io)) ||
            error("$path: unexpected data after the declared k-points")
        return UMatrix(kpoints, data, abspath(path))
    end
end

end
