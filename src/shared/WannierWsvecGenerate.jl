module WannierWsvecGenerate

using LinearAlgebra
using Printf

using ..CrystalCells: wrap_fractional
using ..WannierTypes: HrBlocks, RKey, WsvecEntry, WsvecTable

export generate_wsvec, write_wsvec, validate_wsvec_coverage, assert_wsvec_usable

const WsvecKey = Tuple{RKey, Int, Int}

function _format_key(key::WsvecKey)
    R, i, j = key
    return "($R, $i, $j)"
end

function _sample_join(xs; limit::Int=5)
    isempty(xs) && return ""
    samples = xs[1:min(limit, length(xs))]
    return join(string.(samples), ", ")
end

function _validate_lattice(lattice::AbstractMatrix{<:Real})::Matrix{Float64}
    size(lattice) == (3, 3) || throw(ArgumentError("lattice must be 3 x 3"))
    all(isfinite, lattice) || throw(ArgumentError("lattice contains NaN or Inf"))
    lat = Matrix{Float64}(lattice)
    abs(det(lat)) > eps(Float64) || throw(ArgumentError("lattice must be nonsingular"))
    return lat
end

function _validate_mp_grid(mp_grid::NTuple{3, Int})
    all(>(0), mp_grid) || throw(ArgumentError("mp_grid entries must be positive"))
    return mp_grid
end

function _validate_frozen_axes(frozen_axes::NTuple{3, Bool})
    return frozen_axes
end

function _validate_centers(
    lattice::Matrix{Float64},
    centers_frac::AbstractMatrix{<:Real},
    num_wann::Int,
    centers_cart::Union{Nothing, AbstractMatrix{<:Real}},
    center_consistency_tol::Float64,
)
    size(centers_frac) == (3, num_wann) ||
        throw(ArgumentError("centers_frac must be 3 x num_wann"))
    all(isfinite, centers_frac) || throw(ArgumentError("centers_frac contains NaN or Inf"))

    centers_frac_f = wrap_fractional(Matrix{Float64}(centers_frac))
    centers_cart_from_frac = lattice * centers_frac_f
    if isnothing(centers_cart)
        return centers_frac_f, centers_cart_from_frac
    end

    size(centers_cart) == (3, num_wann) ||
        throw(ArgumentError("centers_cart must be 3 x num_wann"))
    all(isfinite, centers_cart) || throw(ArgumentError("centers_cart contains NaN or Inf"))
    centers_cart_f = Matrix{Float64}(centers_cart)

    maxdiff = maximum(abs.(centers_cart_f .- centers_cart_from_frac))
    if maxdiff > center_consistency_tol
        throw(ArgumentError(
            "centers_cart is inconsistent with lattice * wrap_fractional(centers_frac): " *
            "max difference $maxdiff exceeds tolerance $center_consistency_tol",
        ))
    end
    return centers_frac_f, centers_cart_f
end

function _validate_frozen_axes_compatible(
    hr::HrBlocks,
    frozen_axes::NTuple{3, Bool},
)
    for R in keys(hr.hoppings)
        for dim in 1:3
            if frozen_axes[dim] && R[dim] != 0
                throw(ArgumentError("frozen_axes[$dim] is true but hr contains R=$R"))
            end
        end
    end
    return nothing
end

function _warn_suspicious_R_extent(hr::HrBlocks, mp_grid::NTuple{3, Int})
    samples = RKey[]
    for R in keys(hr.hoppings)
        for dim in 1:3
            limit = ceil(Int, mp_grid[dim] / 2) + 1
            if abs(R[dim]) > limit
                push!(samples, R)
                break
            end
        end
        length(samples) >= 5 && break
    end
    if !isempty(samples)
        @warn (
            "Some hr R-vectors look large for the supplied mp_grid; " *
            "verify that mp_grid matches the Wannier90 calculation."
        ) mp_grid=mp_grid samples=samples
    end
    return nothing
end

function _bvk_translation_candidates(
    mp_grid::NTuple{3, Int},
    search_size::Int,
    frozen_axes::NTuple{3, Bool},
)::Vector{RKey}
    range_for(frozen::Bool) = frozen ? (0:0) : (-search_size:search_size)
    shifts = RKey[]
    for s1 in range_for(frozen_axes[1])
        for s2 in range_for(frozen_axes[2])
            for s3 in range_for(frozen_axes[3])
                push!(shifts, (s1 * mp_grid[1], s2 * mp_grid[2], s3 * mp_grid[3]))
            end
        end
    end
    return shifts
end

function _entry_for_pair(
    R_cart::AbstractVector{<:Real},
    i::Int,
    j::Int,
    centers_cart::AbstractMatrix{<:Real},
    shifts_frac::Vector{RKey},
    shifts_cart::Vector{Vector{Float64}},
    distance_tol::Float64,
)::WsvecEntry
    base_cart = Float64[
        R_cart[1] + centers_cart[1, j] - centers_cart[1, i],
        R_cart[2] + centers_cart[2, j] - centers_cart[2, i],
        R_cart[3] + centers_cart[3, j] - centers_cart[3, i],
    ]

    dmin = Inf
    distances = Vector{Float64}(undef, length(shifts_cart))
    for idx in eachindex(shifts_cart)
        shift = shifts_cart[idx]
        d = sqrt(
            abs2(base_cart[1] + shift[1]) +
            abs2(base_cart[2] + shift[2]) +
            abs2(base_cart[3] + shift[3]),
        )
        distances[idx] = d
        dmin = min(dmin, d)
    end

    tol = max(distance_tol, 16 * eps(max(dmin, 1.0)))
    selected = RKey[]
    for idx in eachindex(distances)
        distances[idx] - dmin <= tol && push!(selected, shifts_frac[idx])
    end
    !isempty(selected) || error("internal error: no wsvec shift selected")

    shift_matrix = Matrix{Int}(undef, 3, length(selected))
    for (idx, shift) in pairs(selected)
        shift_matrix[:, idx] .= shift
    end
    return WsvecEntry(length(selected), shift_matrix)
end

function generate_wsvec(
    hr::HrBlocks,
    lattice::AbstractMatrix{<:Real},
    mp_grid::NTuple{3, Int},
    centers_frac::AbstractMatrix{<:Real};
    centers_cart::Union{Nothing, AbstractMatrix{<:Real}}=nothing,
    search_size::Int=2,
    distance_tol::Float64=1.0e-5,
    frozen_axes::NTuple{3, Bool}=(false, false, false),
    center_consistency_tol::Float64=1.0e-8,
)::WsvecTable
    !isempty(hr.hoppings) || throw(ArgumentError("hr contains no hopping blocks"))
    search_size >= 0 || throw(ArgumentError("search_size must be non-negative"))
    distance_tol >= 0 || throw(ArgumentError("distance_tol must be non-negative"))
    center_consistency_tol >= 0 || throw(ArgumentError("center_consistency_tol must be non-negative"))

    lat = _validate_lattice(lattice)
    _validate_mp_grid(mp_grid)
    _validate_frozen_axes(frozen_axes)
    _, centers_cart_f = _validate_centers(
        lat,
        centers_frac,
        hr.num_wann,
        centers_cart,
        center_consistency_tol,
    )

    _validate_frozen_axes_compatible(hr, frozen_axes)
    _warn_suspicious_R_extent(hr, mp_grid)

    shifts_frac = _bvk_translation_candidates(mp_grid, search_size, frozen_axes)
    shifts_cart = [lat * Float64[s[1], s[2], s[3]] for s in shifts_frac]
    R_cart = Dict{RKey, Vector{Float64}}(
        R => lat * Float64[R[1], R[2], R[3]] for R in keys(hr.hoppings)
    )

    table = Dict{WsvecKey, WsvecEntry}()
    for R in keys(hr.hoppings)
        for j in 1:hr.num_wann, i in 1:hr.num_wann
            table[(R, i, j)] = _entry_for_pair(
                R_cart[R],
                i,
                j,
                centers_cart_f,
                shifts_frac,
                shifts_cart,
                distance_tol,
            )
        end
    end
    return WsvecTable(table)
end

function _duplicate_shifts(entry::WsvecEntry)::Bool
    size(entry.shifts, 2) == entry.n_shift || return false
    seen = Set{RKey}()
    for idx in 1:entry.n_shift
        shift = (entry.shifts[1, idx], entry.shifts[2, idx], entry.shifts[3, idx])
        shift in seen && return true
        push!(seen, shift)
    end
    return false
end

function _negated_shift_set(entry::WsvecEntry)
    shifts = Set{RKey}()
    size(entry.shifts, 2) == entry.n_shift || return shifts
    for idx in 1:entry.n_shift
        push!(shifts, (-entry.shifts[1, idx], -entry.shifts[2, idx], -entry.shifts[3, idx]))
    end
    return shifts
end

function _shift_set(entry::WsvecEntry)
    shifts = Set{RKey}()
    size(entry.shifts, 2) == entry.n_shift || return shifts
    for idx in 1:entry.n_shift
        push!(shifts, (entry.shifts[1, idx], entry.shifts[2, idx], entry.shifts[3, idx]))
    end
    return shifts
end

function validate_wsvec_coverage(
    hr::HrBlocks,
    wsvec::WsvecTable,
)::NamedTuple
    missing_pairs = WsvecKey[]
    bad_entries = Tuple{RKey, Int, Int, String}[]
    extra_pairs = WsvecKey[]
    out_of_bounds_pairs = WsvecKey[]
    duplicate_shifts = WsvecKey[]
    one_sided_pairs = WsvecKey[]
    hermiticity_mismatches = Tuple{RKey, Int, Int, String}[]
    bad_ndegen = RKey[]

    for (R, deg) in hr.ndegen
        deg == 1 || push!(bad_ndegen, R)
    end

    for key in keys(wsvec.table)
        R, i, j = key
        in_bounds = 1 <= i <= hr.num_wann && 1 <= j <= hr.num_wann
        in_bounds || push!(out_of_bounds_pairs, key)
        haskey(hr.hoppings, R) || push!(extra_pairs, key)
    end

    for (R, H) in hr.hoppings
        for j in axes(H, 2), i in axes(H, 1)
            key = (R, i, j)
            entry = get(wsvec.table, key, nothing)
            if isnothing(entry)
                push!(missing_pairs, key)
                continue
            end
            if entry.n_shift < 1 || size(entry.shifts) != (3, entry.n_shift)
                push!(bad_entries, (R, i, j, "invalid shift shape or n_shift"))
                continue
            end
            _duplicate_shifts(entry) && push!(duplicate_shifts, key)
        end
    end

    checked = Set{WsvecKey}()
    for (R, H) in hr.hoppings
        pair_R = (-R[1], -R[2], -R[3])
        for j in axes(H, 2), i in axes(H, 1)
            key = (R, i, j)
            key in checked && continue
            pair_key = (pair_R, j, i)
            push!(checked, key)
            push!(checked, pair_key)

            entry = get(wsvec.table, key, nothing)
            pair_entry = get(wsvec.table, pair_key, nothing)
            if isnothing(entry) || isnothing(pair_entry)
                push!(one_sided_pairs, key)
                continue
            end
            if entry.n_shift < 1 ||
                    pair_entry.n_shift < 1 ||
                    size(entry.shifts) != (3, entry.n_shift) ||
                    size(pair_entry.shifts) != (3, pair_entry.n_shift)
                continue
            end
            _negated_shift_set(entry) == _shift_set(pair_entry) ||
                push!(hermiticity_mismatches, (R, i, j, "shift set is not paired by T -> -T"))
        end
    end

    return (;
        missing_pairs,
        bad_entries,
        extra_pairs,
        out_of_bounds_pairs,
        duplicate_shifts,
        one_sided_pairs,
        hermiticity_mismatches,
        bad_ndegen,
        ok=isempty(missing_pairs) &&
           isempty(bad_entries) &&
           isempty(out_of_bounds_pairs) &&
           isempty(duplicate_shifts),
    )
end

function assert_wsvec_usable(
    hr::HrBlocks,
    wsvec::WsvecTable;
    require_unit_ndegen::Bool=false,
)
    report = validate_wsvec_coverage(hr, wsvec)
    if !isempty(report.missing_pairs)
        error(
            "wsvec is missing $(length(report.missing_pairs)) (R, i, j) " *
            "entr$(length(report.missing_pairs) == 1 ? "y" : "ies") " *
            "(samples: $(_sample_join(_format_key.(report.missing_pairs)))).",
        )
    end
    if !isempty(report.bad_entries)
        error(
            "wsvec has $(length(report.bad_entries)) invalid entr" *
            "$(length(report.bad_entries) == 1 ? "y" : "ies") " *
            "(samples: $(_sample_join(report.bad_entries))).",
        )
    end
    if !isempty(report.out_of_bounds_pairs)
        error(
            "wsvec has $(length(report.out_of_bounds_pairs)) out-of-bounds " *
            "(R, i, j) entr$(length(report.out_of_bounds_pairs) == 1 ? "y" : "ies") " *
            "(samples: $(_sample_join(_format_key.(report.out_of_bounds_pairs)))).",
        )
    end
    if !isempty(report.duplicate_shifts)
        error(
            "wsvec has $(length(report.duplicate_shifts)) entr" *
            "$(length(report.duplicate_shifts) == 1 ? "y" : "ies") with duplicate shifts " *
            "(samples: $(_sample_join(_format_key.(report.duplicate_shifts)))).",
        )
    end
    if require_unit_ndegen && !isempty(report.bad_ndegen)
        error(
            "wsvec requires hr.ndegen == 1 for all R when require_unit_ndegen=true " *
            "(samples: $(_sample_join(report.bad_ndegen))).",
        )
    end
    return report
end

function write_wsvec(
    path::AbstractString,
    wsvec::WsvecTable;
    mp_grid::Union{Nothing, NTuple{3, Int}}=nothing,
    center_policy::Symbol=:unknown,
    header::Union{Nothing, AbstractString}=nothing,
)
    if !isnothing(header) && occursin('\n', String(header))
        throw(ArgumentError("wsvec header must be a single line"))
    end
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        provenance = center_policy == :atomic_assumption ? " atomic-centered-approximation" : ""
        default_header = isnothing(mp_grid) ?
            "## generated by Tools-For-Wannier wsvec$(provenance)" :
            "## generated by Tools-For-Wannier wsvec$(provenance) mp_grid=$(mp_grid)"
        println(io, isnothing(header) ? default_header : String(header))

        for key in sort!(collect(keys(wsvec.table)))
            R, i, j = key
            entry = wsvec.table[key]
            entry.n_shift >= 1 || error("Cannot write wsvec entry $(key): n_shift must be positive")
            size(entry.shifts) == (3, entry.n_shift) ||
                error("Cannot write wsvec entry $(key): shifts must be 3 x n_shift")
            @printf(io, "%5d%5d%5d%5d%5d\n", R[1], R[2], R[3], i, j)
            @printf(io, "%5d\n", entry.n_shift)
            for idx in 1:entry.n_shift
                @printf(io, "%5d%5d%5d\n", entry.shifts[1, idx], entry.shifts[2, idx], entry.shifts[3, idx])
            end
        end
    end
    return path
end

end
