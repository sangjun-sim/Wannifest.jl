module AtomicAngularMomentum

using LinearAlgebra: dot

export P_ORDER, D_ORDER, T2G_ORDER, EG_ORDER
export p_operators, d_operators
export component_expectations

const P_ORDER = ["pz", "px", "py"]
const D_ORDER = ["dz2", "dxz", "dyz", "dx2-y2", "dxy"]
const T2G_ORDER = ["dxz", "dyz", "dxy"]
const EG_ORDER = ["dz2", "dx2-y2"]

function component_expectations(operators, evecs; label::AbstractString, imaginary_tol::Float64=1.0e-10)
    values = zeros(Float64, size(evecs, 2), 5)
    for ib in axes(evecs, 2)
        psi = @view evecs[:, ib]
        x, y, z, squared = map(operators) do op
            value = dot(psi, op * psi)
            abs(imag(value)) <= imaginary_tol ||
                @warn "$label expectation has imaginary residual $(imag(value)); returning real part."
            real(value)
        end
        values[ib, :] .= (x, y, z, sqrt(x^2 + y^2 + z^2), squared)
    end
    return values
end

function p_operators()
    z = 0.0 + 0.0im
    return (
        ComplexF64[z z im; z z z; -im z z],
        ComplexF64[z -im z; im z z; z z z],
        ComplexF64[z z z; z z -im; z im z],
    )
end

function d_operators()
    z = 0.0 + 0.0im
    r3 = sqrt(3.0)
    return (
        ComplexF64[
            z z im*r3 z z
            z z z z im
            -im*r3 z z -im z
            z z im z z
            z -im z z z
        ],
        ComplexF64[
            z -im*r3 z z z
            im*r3 z z -im z
            z z z z -im
            z im z z z
            z z im z z
        ],
        ComplexF64[
            z z z z z
            z z -im z z
            z im z z z
            z z z z -2im
            z z z 2im z
        ],
    )
end

end
