module PairChecks

using ..Hermiticity: assert_hermitian_pair_dict, pair_hermiticity_error

export pair_dict_max_error, check_hr_pair_symmetry

_reverse_key(R::NTuple{3, Int}) = (-R[1], -R[2], -R[3])

function pair_dict_max_error(hops)
    max_err = 0.0
    checked = Set{Tuple{NTuple{3, Int}, NTuple{3, Int}}}()
    for (R, H) in hops
        Rm = _reverse_key(R)
        ordered = isless(R, Rm) ? (R, Rm) : (Rm, R)
        ordered in checked && continue
        push!(checked, ordered)
        Hm = get(hops, Rm, nothing)
        Hm === nothing && return Inf
        max_err = max(max_err, pair_hermiticity_error(H, Hm))
    end
    return max_err
end

check_hr_pair_symmetry(hops; atol::Float64=1e-10) =
    assert_hermitian_pair_dict(hops; atol=atol)

end
