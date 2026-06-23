module AtomicAngularMomentum

export P_ORDER, D_ORDER, T2G_ORDER, EG_ORDER
export p_operators, d_operators, angular_momentum_block

const P_ORDER = ["pz", "px", "py"]
const D_ORDER = ["dz2", "dxz", "dyz", "dx2-y2", "dxy"]
const T2G_ORDER = ["dxz", "dyz", "dxy"]
const EG_ORDER = ["dz2", "dx2-y2"]

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

function _subblock(ops, labels::Vector{String})
    idx = [findfirst(==(label), D_ORDER) for label in labels]
    any(isnothing, idx) && error("internal angular momentum subblock requested non-d orbital labels")
    rows = Int.(idx)
    return (ops[1][rows, rows], ops[2][rows, rows], ops[3][rows, rows])
end

function angular_momentum_block(shell::Symbol, convention::Symbol=:atomic)
    if shell == :p
        convention == :atomic || error("p-shell SOC convention must be atomic")
        Lx, Ly, Lz = p_operators()
        return copy(P_ORDER), Lx, Ly, Lz
    elseif shell == :d
        convention == :atomic || error("d-shell SOC convention must be atomic")
        Lx, Ly, Lz = d_operators()
        return copy(D_ORDER), Lx, Ly, Lz
    elseif shell == :t2g
        convention in (:full_d_projected, :t2g_effective) ||
            error("t2g SOC convention must be full_d_projected or t2g_effective")
        Lx, Ly, Lz = _subblock(d_operators(), T2G_ORDER)
        if convention == :t2g_effective
            Lx, Ly, Lz = -Lx, -Ly, -Lz
        end
        return copy(T2G_ORDER), Lx, Ly, Lz
    end
    error("Unsupported angular momentum shell: $shell")
end

end
