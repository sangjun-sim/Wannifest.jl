module WannierEigIO

export read_eig

"""
    read_eig(path, num_bands, nk) -> Matrix{Float64}

Read a Wannier90 `<seed>.eig` file into a `num_bands × nk` matrix in eV.
Every line holds `band k energy`; bands run fastest and both indices must
count up without gaps, so a missing, repeated or reordered line is an error.
"""
function read_eig(path::AbstractString, num_bands::Integer, nk::Integer)::Matrix{Float64}
    isfile(path) || error("Wannier90 eigenvalue file not found: $(abspath(path))")
    num_bands > 0 && nk > 0 || error("read_eig needs positive num_bands and nk")
    eig = Matrix{Float64}(undef, num_bands, nk)
    expected = 0
    for line in eachline(path)
        fields = split(line)
        isempty(fields) && continue
        expected += 1
        expected <= num_bands * nk ||
            error("$path: more than num_bands*nk = $(num_bands * nk) eigenvalue lines")
        length(fields) == 3 || error("$path: line $expected must hold 'band k energy'")
        band, k = tryparse(Int, fields[1]), tryparse(Int, fields[2])
        energy = tryparse(Float64, replace(fields[3], r"[dD]" => "e"))
        want_k, want_band = fldmod1(expected, Int(num_bands))
        (band, k) == (want_band, want_k) ||
            error("$path: expected band $want_band at k-point $want_k, found '$(fields[1]) $(fields[2])'")
        !isnothing(energy) && isfinite(energy) || error("$path: invalid eigenvalue '$(fields[3])'")
        eig[band, k] = energy
    end
    expected == num_bands * nk ||
        error("$path: found $expected eigenvalue lines, expected num_bands*nk = $(num_bands * nk)")
    return eig
end

end
