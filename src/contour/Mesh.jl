module Mesh

using ..Model: PlaneConfig

export PlaneMesh, generate_plane_mesh, grid_index

struct PlaneMesh
    x_axis::Vector{Float64}
    y_axis::Vector{Float64}
    kpoints::Vector{Vector{Float64}}
    nx::Int
    ny::Int
end

function _axis_values(bounds::Tuple{Float64, Float64}, n::Int)::Vector{Float64}
    n == 1 && return [0.5 * (bounds[1] + bounds[2])]
    return collect(range(bounds[1], bounds[2]; length=n))
end

function grid_index(mesh::PlaneMesh, ik::Integer)
    1 <= ik <= length(mesh.kpoints) || error("mesh index $ik is outside 1:$(length(mesh.kpoints))")
    return (iy=fld(Int(ik) - 1, mesh.nx) + 1, ix=mod1(Int(ik), mesh.nx))
end

function generate_plane_mesh(config::PlaneConfig)::PlaneMesh
    nx, ny = config.mesh
    xs = _axis_values(config.range_x, nx)
    ys = _axis_values(config.range_y, ny)

    kpoints = Vector{Vector{Float64}}(undef, nx * ny)
    ik = 1
    for y in ys, x in xs
        k = zeros(Float64, 3)
        k[config.x_axis] = x
        k[config.y_axis] = y
        k[config.fixed_axis] = config.fixed_value
        kpoints[ik] = k
        ik += 1
    end

    return PlaneMesh(xs, ys, kpoints, nx, ny)
end

end
