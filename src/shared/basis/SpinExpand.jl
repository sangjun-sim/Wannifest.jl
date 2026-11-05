module SpinExpand

using ..SpinLayout
using ..WannierTypes: HrBlocks, RKey

export expand_spinless_hoppings, merge_collinear_hoppings, expand_spinless_basis
export duplicate_hr, merge_collinear_hr

function _output_layout(layout)::Symbol
    mode = SpinLayout.normalize_layout(layout; context="spin expansion layout")
    return mode
end

function _perm_vasp544_to_qe(num_wann::Int)::Vector{Int}
    iseven(num_wann) || error("Expected even num_wann for spinful layout conversion, got $num_wann")
    n = num_wann ÷ 2
    perm = Vector{Int}(undef, num_wann)
    @inbounds for i in 1:n
        perm[2i - 1] = i
        perm[2i] = i + n
    end
    return perm
end

function _from_vasp544_if_needed(H::Matrix{ComplexF64}, layout::Symbol)
    layout == :qe || return H
    perm = _perm_vasp544_to_qe(size(H, 1))
    return H[perm, perm]
end

function _block_diagonal(up::AbstractMatrix{ComplexF64}, dn::AbstractMatrix{ComplexF64})
    nw = size(up, 1)
    size(up) == (nw, nw) || error("up block must be square")
    size(dn) == (nw, nw) || error("down block size mismatch")
    H = zeros(ComplexF64, 2 * nw, 2 * nw)
    H[1:nw, 1:nw] .= up
    H[(nw + 1):(2 * nw), (nw + 1):(2 * nw)] .= dn
    return H
end

function expand_spinless_hoppings(
    hoppings::Dict{RKey, Matrix{ComplexF64}},
    num_wann::Integer,
    layout,
)::Dict{RKey, Matrix{ComplexF64}}
    nw = Int(num_wann)
    mode = _output_layout(layout)
    out = Dict{RKey, Matrix{ComplexF64}}()
    for (R, H0) in hoppings
        size(H0) == (nw, nw) || error("Block for R=$R does not have size ($nw, $nw)")
        out[R] = _from_vasp544_if_needed(_block_diagonal(H0, H0), mode)
    end
    return out
end

function merge_collinear_hoppings(
    up_hoppings::Dict{RKey, Matrix{ComplexF64}},
    dn_hoppings::Dict{RKey, Matrix{ComplexF64}},
    num_wann::Integer,
    layout,
)::Dict{RKey, Matrix{ComplexF64}}
    nw = Int(num_wann)
    keys_up = sort!(collect(keys(up_hoppings)))
    keys_dn = sort!(collect(keys(dn_hoppings)))
    keys_up == keys_dn || error("R-point set mismatch between up/down hr.dat inputs")
    mode = _output_layout(layout)
    out = Dict{RKey, Matrix{ComplexF64}}()
    for R in keys_up
        Hup = up_hoppings[R]
        Hdn = dn_hoppings[R]
        size(Hup) == (nw, nw) || error("Up block for R=$R does not have size ($nw, $nw)")
        size(Hdn) == (nw, nw) || error("Down block for R=$R does not have size ($nw, $nw)")
        out[R] = _from_vasp544_if_needed(_block_diagonal(Hup, Hdn), mode)
    end
    return out
end

function expand_spinless_basis(entries, layout; entry_constructor)
    mode = _output_layout(layout)
    out = Any[]
    if mode == :vasp544
        for spin in (:up, :dn), entry in entries
            push!(out, entry_constructor(length(out) + 1, entry, spin))
        end
    else
        for entry in entries, spin in (:up, :dn)
            push!(out, entry_constructor(length(out) + 1, entry, spin))
        end
    end
    return out
end

function duplicate_hr(hr::HrBlocks, layout)::HrBlocks
    hops = expand_spinless_hoppings(hr.hoppings, hr.num_wann, layout)
    return HrBlocks(hr.header, 2 * hr.num_wann, hr.nrpts, hops, copy(hr.ndegen), hr.normalization)
end

function merge_collinear_hr(hr_up::HrBlocks, hr_dn::HrBlocks, layout)::HrBlocks
    hr_up.num_wann == hr_dn.num_wann ||
        error("num_wann mismatch: $(hr_up.num_wann) vs $(hr_dn.num_wann)")
    hr_up.nrpts == hr_dn.nrpts ||
        error("nrpts mismatch: $(hr_up.nrpts) vs $(hr_dn.nrpts)")
    for R in keys(hr_up.hoppings)
        hr_up.ndegen[R] == hr_dn.ndegen[R] ||
            error("ndegen mismatch at R=$R: $(hr_up.ndegen[R]) vs $(hr_dn.ndegen[R])")
    end
    hops = merge_collinear_hoppings(hr_up.hoppings, hr_dn.hoppings, hr_up.num_wann, layout)
    header = hr_up.header == hr_dn.header ? hr_up.header :
        string(hr_up.header, " | merged with | ", hr_dn.header)
    return HrBlocks(header, 2 * hr_up.num_wann, hr_up.nrpts, hops, copy(hr_up.ndegen), hr_up.normalization)
end

end
