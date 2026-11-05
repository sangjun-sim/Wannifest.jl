module BandWindow

using LinearAlgebra: Diagonal

using ..WannierUmatIO: UMatrix

export kept_bands, outer_window, aligned_u, subspace_offdiagonal
export window_edge_degeneracies, check_eig_against_qe

const ZERO_ROW_TOL = 1e-8

"""QE band index of every band Wannier90 kept (`exclude_bands` removed, length `num_bands`)."""
function kept_bands(nbnd_qe::Integer, exclude::Vector{Int}, num_bands::Integer)::Vector{Int}
    all(band -> 1 <= band <= nbnd_qe, exclude) ||
        error("exclude_bands holds an index outside the QE band range 1:$nbnd_qe")
    kept = setdiff(1:nbnd_qe, exclude)
    length(kept) == num_bands ||
        error("QE has $nbnd_qe bands and exclude_bands removes $(length(exclude)), leaving $(length(kept)); Wannier90 num_bands is $num_bands")
    return kept
end

"""
    outer_window(eig, win_min, win_max, num_wann) -> (nfirstwin, ndimwin)

Rebuild the outer window of every k-point with the rule of Wannier90 `dis_windows`:
the window starts at the first band with `win_min ≤ E ≤ win_max` and ends at the
last band with `E ≤ win_max`. `eig` is `num_bands × nk`, ascending in each column.
"""
function outer_window(eig::Matrix{Float64}, win_min::Real, win_max::Real, num_wann::Integer)
    win_min <= win_max || error("dis_win_min ($win_min) exceeds dis_win_max ($win_max)")
    nk = size(eig, 2)
    nfirstwin, ndimwin = Vector{Int}(undef, nk), Vector{Int}(undef, nk)
    for k in 1:nk
        energies = view(eig, :, k)
        issorted(energies) || error("Eigenvalues at k-point $k are not in ascending order")
        imin = findfirst(e -> win_min <= e <= win_max, energies)
        isnothing(imin) &&
            error("The outer energy window [$win_min, $win_max] eV contains no eigenvalues at k-point $k")
        imax = findlast(<=(win_max), energies)
        nfirstwin[k], ndimwin[k] = imin, imax - imin + 1
        ndimwin[k] >= num_wann ||
            error("The outer window holds $(ndimwin[k]) states at k-point $k, fewer than num_wann = $num_wann")
    end
    return nfirstwin, ndimwin
end

"""
    aligned_u(u_dis, u, nfirstwin, ndimwin) -> Array{ComplexF64, 3}

Ũ(k) = U_dis(k)·U(k) as `num_bands × num_wann`, with every row at its band index.
Row m of `_u_dis.mat` is the m-th band *inside the outer window*, i.e. band
`nfirstwin(k)+m−1`, and Wannier90 zeroes the rows past `ndimwin(k)`. A nonzero
row there means the rebuilt window differs from the one Wannier90 used.
"""
function aligned_u(u_dis::UMatrix, u::UMatrix, nfirstwin::Vector{Int}, ndimwin::Vector{Int})
    nb, nw, nk = size(u_dis.data)
    size(u.data) == (nw, nw, nk) ||
        error("$(u.source) must hold $nw×$nw matrices at $nk k-points to match $(u_dis.source)")
    length(nfirstwin) == nk && length(ndimwin) == nk || error("Outer window arrays must have one entry per k-point")
    aligned = zeros(ComplexF64, nb, nw, nk)
    for k in 1:nk
        d = ndimwin[k]
        rows = nfirstwin[k]:(nfirstwin[k] + d - 1)
        1 <= first(rows) && last(rows) <= nb || error("Outer window $rows at k-point $k leaves the band range 1:$nb")
        maximum(abs, @view u_dis.data[(d + 1):end, :, k]; init=0.0) <= ZERO_ROW_TOL ||
            error("U_dis rows outside the reconstructed outer window at k=$k in $(u_dis.source); " *
                  "the .eig file or the dis_win_min/max values differ from the Wannier90 run")
        aligned[rows, :, k] = u_dis.data[1:d, :, k] * u.data[:, :, k]
    end
    return aligned
end

"""
    subspace_offdiagonal(u_dis, eig, nfirstwin, ndimwin) -> (value, k)

Largest off-diagonal element (eV) of `U_dis† diag(E_window) U_dis` and its k-point.
Wannier90 rotates `u_matrix_opt` until this matrix is diagonal, so a large value
flags a window shifted against the Wannier90 run, which the zero-row test of
`aligned_u` cannot see when the rebuilt window is too wide.
"""
function subspace_offdiagonal(u_dis::UMatrix, eig::Matrix{Float64}, nfirstwin::Vector{Int}, ndimwin::Vector{Int})
    worst, worst_k = 0.0, 0
    for k in axes(u_dis.data, 3)
        d = ndimwin[k]
        block = @view u_dis.data[1:d, :, k]
        energies = @view eig[nfirstwin[k]:(nfirstwin[k] + d - 1), k]
        h = block' * Diagonal(energies) * block
        offdiagonal = maximum(abs(h[i, j]) for i in axes(h, 1), j in axes(h, 2) if i != j; init=0.0)
        offdiagonal > worst && ((worst, worst_k) = (offdiagonal, k))
    end
    return worst, worst_k
end

"""
k-points where an edge of the outer window cuts a (near-)degenerate multiplet:
`|E(imin−1) − E(imin)|` or `|E(imax) − E(imax+1)|` below `degeneracy_tol` (eV).
The subspace itself is then ambiguous, so T is not protected by gauge invariance.
"""
function window_edge_degeneracies(eig::Matrix{Float64}, nfirstwin::Vector{Int}, ndimwin::Vector{Int}; degeneracy_tol::Real)
    flagged = Int[]
    for k in axes(eig, 2)
        imin, imax = nfirstwin[k], nfirstwin[k] + ndimwin[k] - 1
        lower = imin > 1 && abs(eig[imin, k] - eig[imin - 1, k]) < degeneracy_tol
        upper = imax < size(eig, 1) && abs(eig[imax + 1, k] - eig[imax, k]) < degeneracy_tol
        (lower || upper) && push!(flagged, k)
    end
    return flagged
end

"""
    check_eig_against_qe(eig, qe_eig, kept, kindex; eig_tol, label) -> max|ΔE|

Compare `<seed>.eig` (`num_bands × nk`, eV) with the QE XML eigenvalues of the kept
bands; `kindex[k]` is the QE k index of Wannier k-point `k`. A difference above
`eig_tol` reports the offending (k, band, ΔE): it means another nscf outdir, a wrong
`exclude_bands`, a wrong k correspondence or swapped LSDA channels.
"""
function check_eig_against_qe(eig::Matrix{Float64}, qe_eig::Matrix{Float64}, kept::Vector{Int}, kindex::Vector{Int};
                              eig_tol::Real, label::AbstractString)
    size(eig) == (length(kept), length(kindex)) || error("$label: .eig must be num_bands × nk")
    worst, worst_k, worst_band = 0.0, 0, 0
    for k in axes(eig, 2), (band, qe_band) in enumerate(kept)
        delta = abs(eig[band, k] - qe_eig[qe_band, kindex[k]])
        delta > worst && ((worst, worst_k, worst_band) = (delta, k, band))
    end
    worst <= eig_tol ||
        error("$label: .eig differs from the QE XML eigenvalues at k-point $worst_k, band $worst_band " *
              "(QE band $(kept[worst_band]), QE k-point $(kindex[worst_k])): |ΔE| = $worst eV > eig_tol = $eig_tol eV. " *
              "Check the nscf outdir, exclude_bands, the k-point order and the up/down seeds.")
    return worst
end

end
