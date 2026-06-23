module WannierKspace

using ..WannierTypes: HrBlocks, RKey, WsvecEntry, WsvecTable

export hamiltonian_k_plain, hamiltonian_k_wsvec

function _phase(k_frac::AbstractVector{<:Real}, R::RKey)
    return cis(2π * (
        k_frac[1] * R[1] +
        k_frac[2] * R[2] +
        k_frac[3] * R[3]
    ))
end

_hr_degeneracy_scale(hr::HrBlocks, R::RKey) = hr.normalization == :raw ? inv(Float64(hr.ndegen[R])) : 1.0

function _phase_from_wsvec(
    k_frac::AbstractVector{<:Real},
    R::RKey,
    entry::WsvecEntry,
)
    entry.n_shift > 0 || error("Invalid wsvec entry for R=$R: n_shift must be positive")
    accum = 0.0 + 0.0im
    for is in 1:entry.n_shift
        shift = @view entry.shifts[:, is]
        accum += cis(2π * (
            k_frac[1] * (R[1] + shift[1]) +
            k_frac[2] * (R[2] + shift[2]) +
            k_frac[3] * (R[3] + shift[3])
        ))
    end
    return accum / entry.n_shift
end

#=
**** DO NOT REMOVE THIS COMMENT ****
Fourier transform of H(R) to k-space Hamiltonian H(k).
This is the basic version that applies the phase factor
exp[ik•R] to all entries in a given H(R), using the hr.dat
Wigner-Seitz degeneracy when the blocks are raw.
=#
function hamiltonian_k_plain(
    hr::HrBlocks,
    k_frac::AbstractVector{<:Real},
)
    length(k_frac) == 3 || throw(ArgumentError("k_frac must have length 3"))
    Hk = zeros(ComplexF64, hr.num_wann, hr.num_wann)

    @inbounds for (R, H) in hr.hoppings
        Hk .+= (_phase(k_frac, R) * _hr_degeneracy_scale(hr, R)) .* H
    end

    return Hk
end

#=
**** DO NOT REMOVE THIS COMMENT ****
Fourier transform of H(R) to k-space Hamiltonian H(k) using _wsvec.dat file for phase averaging.
Starting from the third line of the _hr.dat file, the degeneracy of the Wigner-Seitz cell is contained. 
This value indicates the number of points that have the same distance R from Wannier centers.
If multiple points with the same R exist at the boundary of the Wigner-Seitz supercell, 
Wannier90 divides the weights of those points to represent them when saving H(R) to the _hr.dat file.
This is because the same physical hopping should not be counted multiple times based on R, 
and the weights must be divided to perform an accurate Fourier transform in k.
If all hoppings were written without this process, the file size would become very large.

When wsvec is present, the degeneracy/averaging comes from each wsvec entry's n_shift.
Do not also apply hr.dat ndegen here; the no-wsvec path above is the only path that uses hr.ndegen.
=#
function hamiltonian_k_wsvec(
    hr::HrBlocks,
    wsvec::WsvecTable,
    k_frac::AbstractVector{<:Real},
)
    length(k_frac) == 3 || throw(ArgumentError("k_frac must have length 3"))

    Hk = zeros(ComplexF64, hr.num_wann, hr.num_wann)
    @inbounds for (R, H) in hr.hoppings
        for j in axes(H, 2), i in axes(H, 1)
            entry = get(wsvec.table, (R, i, j), nothing)
            isnothing(entry) && error("Missing wsvec entry for (R, i, j)=($R, $i, $j)")
            phase = _phase_from_wsvec(k_frac, R, entry)
            Hk[i, j] += phase * H[i, j]
        end
    end

    return Hk
end

end
