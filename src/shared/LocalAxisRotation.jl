module LocalAxisRotation

using LinearAlgebra

using ..BasisLabelNormalize: canonical_orbital, canonical_spin

export LocalBasisEntry, AxisSpec, parse_axes
export basis_entries_from_win, build_rotation_transform

const P_ORDER = ["pz", "px", "py"]
const D_ORDER = ["dz2", "dxz", "dyz", "dx2-y2", "dxy"]
const T2G_ORDER = ["dxz", "dyz", "dxy"]
const D_NORM2 = 1.5

struct LocalBasisEntry
    index::Int
    site::String
    orbital::String
    spin::String
end

struct AxisSpec
    site::String
    source_x::Vector{Float64}
    source_z::Vector{Float64}
    target_x::Vector{Float64}
    target_z::Vector{Float64}
end

function basis_entries_from_win(basis)
    return LocalBasisEntry[
        LocalBasisEntry(orb.index, orb.site_label, canonical_orbital(orb.orbital), canonical_spin(orb.spin))
        for orb in basis.orbitals
    ]
end

function _vec3(raw, context::AbstractString)
    raw isa AbstractVector && length(raw) == 3 || error("$context must be a length-3 array")
    return [Float64(raw[1]), Float64(raw[2]), Float64(raw[3])]
end

function _parse_axis_record(raw, index::Integer)
    raw isa AbstractVector ||
        error("local_axes entry $index must be [site, source_z, source_x, target_z, target_x]")
    length(raw) == 5 ||
        error("local_axes entry $index must have 5 items: [site, source_z, source_x, target_z, target_x]")
    raw[1] isa AbstractString || error("local_axes entry $index item 1 must be a site string")
    site = strip(String(raw[1]))
    isempty(site) && error("local_axes entry $index site cannot be empty")
    source_z = _vec3(raw[2], "source_z for site $site")
    source_x = _vec3(raw[3], "source_x for site $site")
    target_z = _vec3(raw[4], "target_z for site $site")
    target_x = _vec3(raw[5], "target_x for site $site")
    return AxisSpec(site, source_x, source_z, target_x, target_z)
end

function parse_axes(raw_axes)
    raw_axes isa AbstractVector ||
        error("local_axes must be an array of [site, source_z, source_x, target_z, target_x] records")
    axes = AxisSpec[_parse_axis_record(raw, i) for (i, raw) in enumerate(raw_axes)]
    isempty(axes) && error("projection.basis_rotation.local_axes has no entries")
    length(unique(axis.site for axis in axes)) == length(axes) ||
        error("projection.basis_rotation.local_axes has duplicate site entries")
    return axes
end

function _normalize_axis(v::AbstractVector{<:Real}, context::AbstractString)
    n = norm(v)
    n > 0 || error("$context must be nonzero")
    return Float64.(v) ./ n
end

function local_frame(xaxis::AbstractVector{<:Real}, zaxis::AbstractVector{<:Real}; context::AbstractString="")
    zhat = _normalize_axis(zaxis, "$context z-axis")
    xraw = _normalize_axis(xaxis, "$context x-axis")
    dot_xz = dot(xraw, zhat)
    xproj = xraw .- dot_xz .* zhat
    norm(xproj) > 1.0e-12 || error("$context x-axis is parallel to z-axis")
    abs(dot_xz) <= 1.0e-8 || @warn "$context x-axis is not orthogonal to z-axis; using its perpendicular projection"
    xhat = xproj ./ norm(xproj)
    yhat = cross(zhat, xhat)
    return hcat(xhat, yhat, zhat)
end

function _p_axis(label::AbstractString)
    label == "px" && return 1
    label == "py" && return 2
    label == "pz" && return 3
    error("Unsupported p orbital label: $label")
end

function p_rotation(source_frame::AbstractMatrix{<:Real}, target_frame::AbstractMatrix{<:Real})
    U = zeros(Float64, length(P_ORDER), length(P_ORDER))
    for (a, target_label) in enumerate(P_ORDER)
        target_vec = target_frame[:, _p_axis(target_label)]
        for (i, source_label) in enumerate(P_ORDER)
            source_vec = source_frame[:, _p_axis(source_label)]
            U[a, i] = dot(target_vec, source_vec)
        end
    end
    return U
end

function _d_local_mats()
    rt3 = sqrt(3.0)
    return Dict{String, Matrix{Float64}}(
        "dz2" => [-0.5 0.0 0.0; 0.0 -0.5 0.0; 0.0 0.0 1.0],
        "dxz" => [0.0 0.0 rt3 / 2; 0.0 0.0 0.0; rt3 / 2 0.0 0.0],
        "dyz" => [0.0 0.0 0.0; 0.0 0.0 rt3 / 2; 0.0 rt3 / 2 0.0],
        "dx2-y2" => [rt3 / 2 0.0 0.0; 0.0 -rt3 / 2 0.0; 0.0 0.0 0.0],
        "dxy" => [0.0 rt3 / 2 0.0; rt3 / 2 0.0 0.0; 0.0 0.0 0.0],
    )
end

const D_LOCAL_MATS = _d_local_mats()

function d_rotation(source_frame::AbstractMatrix{<:Real}, target_frame::AbstractMatrix{<:Real})
    source_global = Dict(label => source_frame * D_LOCAL_MATS[label] * transpose(source_frame) for label in D_ORDER)
    target_global = Dict(label => target_frame * D_LOCAL_MATS[label] * transpose(target_frame) for label in D_ORDER)
    U = zeros(Float64, length(D_ORDER), length(D_ORDER))
    for (a, target_label) in enumerate(D_ORDER)
        for (i, source_label) in enumerate(D_ORDER)
            U[a, i] = sum(target_global[target_label] .* source_global[source_label]) / D_NORM2
        end
    end
    return U
end

function _entries_by_site(entries::Vector{LocalBasisEntry})
    by_site = Dict{String, Vector{LocalBasisEntry}}()
    for entry in entries
        push!(get!(by_site, entry.site, LocalBasisEntry[]), entry)
    end
    return by_site
end

function _entries_by_spin(entries::Vector{LocalBasisEntry})
    by_spin = Dict{String, Vector{LocalBasisEntry}}()
    for entry in entries
        push!(get!(by_spin, entry.spin, LocalBasisEntry[]), entry)
    end
    return by_spin
end

function _entries_by_orbital(entries::Vector{LocalBasisEntry}, site::AbstractString, spin::AbstractString)
    by_orbital = Dict{String, LocalBasisEntry}()
    for entry in entries
        haskey(by_orbital, entry.orbital) &&
            error("site $site spin $spin has duplicate orbital $(entry.orbital)")
        by_orbital[entry.orbital] = entry
    end
    return by_orbital
end

_has_all(by_orbital, labels::Vector{String}) = all(label -> haskey(by_orbital, label), labels)
_has_any(by_orbital, labels::Vector{String}) = any(label -> haskey(by_orbital, label), labels)

function _place_block!(U::Matrix{ComplexF64}, labels::Vector{String}, by_orbital, block::AbstractMatrix{<:Real})
    indices = [by_orbital[label].index for label in labels]
    U[indices, indices] .= ComplexF64.(block)
    return nothing
end

function _d_index(label::AbstractString)
    idx = findfirst(==(String(label)), D_ORDER)
    isnothing(idx) && error("Unsupported d orbital label: $label")
    return idx
end

function _apply_site_spin_rotation!(
    U::Matrix{ComplexF64},
    entries::Vector{LocalBasisEntry},
    axis::AxisSpec,
    source_frame::Matrix{Float64},
    target_frame::Matrix{Float64},
    strict_t2g::Bool,
    leakage_tol::Float64,
)
    by_orbital = _entries_by_orbital(entries, axis.site, first(entries).spin)

    if _has_all(by_orbital, P_ORDER)
        _place_block!(U, P_ORDER, by_orbital, p_rotation(source_frame, target_frame))
    elseif _has_any(by_orbital, P_ORDER)
        missing = setdiff(P_ORDER, collect(keys(by_orbital)))
        error("site $(axis.site) spin $(first(entries).spin) has an incomplete p shell; missing $(join(missing, ", "))")
    end

    d_present = [label for label in D_ORDER if haskey(by_orbital, label)]
    if _has_all(by_orbital, D_ORDER)
        _place_block!(U, D_ORDER, by_orbital, d_rotation(source_frame, target_frame))
    elseif Set(d_present) == Set(T2G_ORDER)
        full = d_rotation(source_frame, target_frame)
        t2g_rows = [_d_index(label) for label in T2G_ORDER]
        eg_rows = [_d_index(label) for label in ("dz2", "dx2-y2")]
        leakage = maximum(sqrt(sum(abs2, full[eg_rows, col])) for col in t2g_rows)
        if leakage > leakage_tol
            msg = "site $(axis.site) spin $(first(entries).spin) t2g rotation leaks into eg by $(leakage); projecting back to t2g-only block"
            strict_t2g ? error(msg) : @warn msg
        end
        _place_block!(U, T2G_ORDER, by_orbital, full[t2g_rows, t2g_rows])
    elseif !isempty(d_present)
        missing_t2g = setdiff(T2G_ORDER, collect(keys(by_orbital)))
        missing_d = setdiff(D_ORDER, collect(keys(by_orbital)))
        error(
            "site $(axis.site) spin $(first(entries).spin) has unsupported partial d shell; " *
            "missing t2g=$(join(missing_t2g, ", ")), missing full d=$(join(missing_d, ", "))"
        )
    end

    return nothing
end

function build_rotation_transform(
    num_wann::Integer,
    basis::Vector{LocalBasisEntry},
    axes::Vector{AxisSpec};
    strict_t2g::Bool=false,
    leakage_tol::Float64=1.0e-8,
)
    nw = Int(num_wann)
    all(entry -> 1 <= entry.index <= nw, basis) ||
        error("basis index outside hr range 1:$nw")

    U = Matrix{ComplexF64}(I, nw, nw)
    by_site = _entries_by_site(basis)
    for axis in axes
        site_entries = get(by_site, axis.site, nothing)
        isnothing(site_entries) && error("local_axes references site $(axis.site), but basis has no matching entries")
        source_frame = local_frame(axis.source_x, axis.source_z; context="site $(axis.site) source")
        target_frame = local_frame(axis.target_x, axis.target_z; context="site $(axis.site) target")
        for (_, spin_entries) in sort!(collect(_entries_by_spin(site_entries)); by=first)
            _apply_site_spin_rotation!(
                U,
                spin_entries,
                axis,
                source_frame,
                target_frame,
                strict_t2g,
                leakage_tol,
            )
        end
    end
    return U
end

end
