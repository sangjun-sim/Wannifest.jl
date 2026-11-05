module HrHermiticity

using ..WannierTypes: RKey

export complete_hermiticity!, reverse_key

reverse_key(R::RKey)::RKey = (-R[1], -R[2], -R[3])

function _complete_zero_block!(H::Matrix{ComplexF64}; atol::Float64=1e-12)
    n = size(H, 1)
    for i in 1:n
        abs(imag(H[i, i])) <= atol ||
            error("Onsite diagonal element at ($i,$i) has non-zero imaginary part: $(H[i, i])")
        H[i, i] = ComplexF64(real(H[i, i]), 0.0)
        for j in (i + 1):n
            a = H[i, j]
            b = H[j, i]
            a_zero = abs(a) <= atol
            b_zero = abs(b) <= atol
            if a_zero && b_zero
                continue
            elseif a_zero
                H[i, j] = conj(b)
            elseif b_zero
                H[j, i] = conj(a)
            elseif !isapprox(a, conj(b); atol=atol, rtol=0.0)
                error("R=(0,0,0) block is not Hermitian at ($i,$j): $a vs $(conj(b))")
            end
        end
    end
    H .= 0.5 .* (H .+ H')
    return H
end

function _complete_partner_blocks!(H::Matrix{ComplexF64}, Hm::Matrix{ComplexF64}; atol::Float64=1e-12)
    n = size(H, 1)
    size(Hm) == (n, n) || error("Hermitian partner block dimension mismatch")
    for i in 1:n, j in 1:n
        a = H[i, j]
        b = conj(Hm[j, i])
        a_zero = abs(a) <= atol
        b_zero = abs(b) <= atol
        if a_zero && b_zero
            continue
        elseif a_zero
            H[i, j] = b
        elseif b_zero
            Hm[j, i] = conj(a)
        elseif !isapprox(a, b; atol=atol, rtol=0.0)
            error("Hermitian partner mismatch at ($i,$j): $a vs $b")
        end
        avg = 0.5 * (H[i, j] + conj(Hm[j, i]))
        H[i, j] = avg
        Hm[j, i] = conj(avg)
    end
    return nothing
end

function complete_hermiticity!(hops::Dict{RKey, Matrix{ComplexF64}}; atol::Float64=1e-12)
    for (R, H) in collect(hops)
        if R == (0, 0, 0)
            _complete_zero_block!(H; atol=atol)
            continue
        end
        Rm = reverse_key(R)
        if !haskey(hops, Rm)
            hops[Rm] = copy(H')
            continue
        end
        _complete_partner_blocks!(H, hops[Rm]; atol=atol)
    end
    return hops
end

end
