module Hermiticity

export pair_hermiticity_error, assert_hermitian_pair_dict

function pair_hermiticity_error(A::AbstractMatrix{<:Number}, B::AbstractMatrix{<:Number})
    return maximum(abs, A .- B')
end

function assert_hermitian_pair_dict(hops; atol::Float64=1e-10)
    max_err = 0.0
    for (R, H) in hops
        Rm = (-R[1], -R[2], -R[3])
        Hm = get(hops, Rm, nothing)
        Hm === nothing && error("Missing Hermitian partner hopping for R=$R")
        err = pair_hermiticity_error(H, Hm)
        max_err = max(max_err, err)
    end
    max_err <= atol || error("Hermiticity check failed with max error $max_err > $atol")
    return max_err
end

end
