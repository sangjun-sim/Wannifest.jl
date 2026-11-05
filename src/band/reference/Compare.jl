module Compare

using LinearAlgebra
using Printf

using ..KPath: KPathResult, format_point_label
using ..Model: BandComparison, BandResult, ReferenceConfig
using ..ReferenceData: ReferenceBandData

export qe_kpath, aligned_energies, reference_fermi_level, compare_bands, combined_rms, write_comparison

const KPOINT_TOL = 1e-10

_same_kpoint(left, right) = norm(left .- right) <= KPOINT_TOL

function _is_corner(left, right)
    left_norm = norm(left)
    right_norm = norm(right)
    (left_norm <= KPOINT_TOL || right_norm <= KPOINT_TOL) && return false
    return abs(dot(left, right) / (left_norm * right_norm) - 1.0) > 1e-8
end

function _segments(kpoints::Vector{Vector{Float64}})
    ranges = UnitRange{Int}[]
    start = 1
    index = 2
    while index < length(kpoints)
        if _same_kpoint(kpoints[index], kpoints[index + 1])
            push!(ranges, start:index)
            start = index + 1
            index += 1
        elseif _is_corner(kpoints[index] .- kpoints[index - 1], kpoints[index + 1] .- kpoints[index])
            push!(ranges, start:index)
            start = index
        end
        index += 1
    end
    push!(ranges, start:length(kpoints))
    return ranges
end

function _tick_indices(ranges::Vector{UnitRange{Int}})
    indices = Int[first(ranges[1])]
    for range in ranges
        first(range) != indices[end] && push!(indices, first(range))
        push!(indices, last(range))
    end
    return indices
end

const PATH_BREAK = "|"

_is_break_label(label::AbstractString) = occursin(PATH_BREAK, label)

"""
    _tick_plan(tick_indices, labels)

Match `reference.k_labels` against the path endpoints `_segments` found in the
QE k-point list. A label written with `|` (`"B_1|B"`) marks a discontinuity, so
it stands for two endpoints: the last k-point of the line that ends and the
first k-point of the line that starts. Returns the k-point index each label
ticks at, plus the indices that open a line across such a jump.
"""
function _tick_plan(tick_indices::Vector{Int}, labels::Vector{String})
    breaks = count(_is_break_label, labels)
    expected = length(labels) + breaks
    expected == length(tick_indices) || error(
        "reference.k_labels describes $expected path endpoints " *
        "($(length(labels)) labels, $breaks discontinuity markers), " *
        "but QE path has $(length(tick_indices)) endpoints. " *
        "Write a discontinuity as one label containing \"$PATH_BREAK\", e.g. \"B_1|B\"",
    )

    tick_at = Int[]
    jump_starts = Set{Int}()
    cursor = 1
    for (ilabel, label) in enumerate(labels)
        index = tick_indices[cursor]
        cursor += 1
        if _is_break_label(label)
            next_index = tick_indices[cursor]
            cursor += 1
            next_index == index + 1 || error(
                "reference.k_labels[$ilabel] = \"$label\" marks a discontinuity, but the QE " *
                "k-points it falls between ($index and $next_index) are not adjacent",
            )
            push!(jump_starts, next_index)
        end
        push!(tick_at, index)
    end
    return tick_at, jump_starts
end

# The straight hop across a discontinuity is not a band segment, so it must not
# be drawn. `_segments` reports it as the two-point range that ends where the
# next line opens; a repeated k-point leaves no such range to drop.
_is_jump_range(range::UnitRange{Int}, jump_starts::Set{Int}) =
    last(range) in jump_starts && first(range) == last(range) - 1

function qe_kpath(qe::ReferenceBandData, lattice, labels::Vector{String})::KPathResult
    kpoints = qe.kpoints_frac
    length(kpoints) >= 2 || error("QE reference requires at least two k-points")
    ranges = _segments(kpoints)
    tick_at, jump_starts = _tick_plan(_tick_indices(ranges), labels)

    for index in 2:length(ranges)
        first(ranges[index]) > last(ranges[index - 1]) &&
            push!(jump_starts, first(ranges[index]))
    end
    distances = zeros(Float64, length(kpoints))
    for index in 2:length(kpoints)
        step = index in jump_starts ? 0.0 :
            norm(lattice.reciprocal_lattice * (kpoints[index] .- kpoints[index - 1]))
        distances[index] = distances[index - 1] + step
    end
    plot_ranges = filter(range -> !_is_jump_range(range, jump_starts), ranges)
    formatted_labels = format_point_label.(labels)
    return KPathResult(kpoints, distances, distances[tick_at], formatted_labels, plot_ranges, true)
end

function _vbm(eigenvalues::AbstractArray{<:Real}, fermi_energy::Real)
    occupied = eigenvalues[eigenvalues .<= fermi_energy]
    isempty(occupied) && error("cannot align to VBM: no eigenvalue is <= Fermi energy")
    return maximum(occupied)
end

"""
    reference_fermi_level(qe, reference)

Reference level expressed in the aligned energies that `aligned_energies`
returns, so a comparison plot can mark it with the dashed line. It is the
plotted zero for `fermi` and `vbm` alignment (the Fermi level and the VBM
respectively) and the raw Fermi energy for `none`.
"""
function reference_fermi_level(qe::ReferenceBandData, reference::ReferenceConfig)::Float64
    reference.align == :fermi && return 0.0
    reference.align == :vbm && return 0.0
    reference.align == :none && return qe.fermi_energy
    error("unsupported reference.align: $(reference.align)")
end

function aligned_energies(
    wannier::BandResult,
    qe::ReferenceBandData,
    reference::ReferenceConfig,
)
    wannier_values = copy(wannier.eigenvalues)
    qe_values = copy(qe.eigenvalues)
    if reference.align == :fermi
        wannier_values .-= qe.fermi_energy
        qe_values .-= qe.fermi_energy
    elseif reference.align == :vbm
        wannier_values .-= _vbm(wannier_values, qe.fermi_energy)
        qe_values .-= _vbm(qe_values, qe.fermi_energy)
    end
    return wannier_values, qe_values
end

function aligned_energies(
    up::BandResult,
    down::Union{Nothing, BandResult},
    qe::ReferenceBandData,
    reference::ReferenceConfig,
)
    isnothing(down) && return aligned_energies(up, qe, reference)
    qe_down = qe.eigenvalues_down
    isnothing(qe_down) && error("QE reference is not lsda=true, but two Wannier channels were provided")

    up_values = copy(up.eigenvalues)
    down_values = copy(down.eigenvalues)
    qe_up_values = copy(qe.eigenvalues)
    qe_down_values = copy(qe_down)
    if reference.align == :fermi
        up_values .-= qe.fermi_energy
        down_values .-= qe.fermi_energy
        qe_up_values .-= qe.fermi_energy
        qe_down_values .-= qe.fermi_energy
    elseif reference.align == :vbm
        wannier_vbm = _vbm(vcat(vec(up_values), vec(down_values)), qe.fermi_energy)
        qe_vbm = _vbm(vcat(vec(qe_up_values), vec(qe_down_values)), qe.fermi_energy)
        up_values .-= wannier_vbm
        down_values .-= wannier_vbm
        qe_up_values .-= qe_vbm
        qe_down_values .-= qe_vbm
    end
    return up_values, down_values, qe_up_values, qe_down_values
end

function compare_bands(
    wannier::BandResult,
    qe_values::AbstractMatrix{<:Real},
    window::Tuple{Float64, Float64},
)::BandComparison
    nkpoints, nwannier = size(wannier.eigenvalues)
    size(qe_values, 1) == nkpoints || error("QE and Wannier k-point counts differ")

    count_in_window = zeros(Int, nwannier)
    rms = fill(NaN, nwannier)
    max_error = fill(NaN, nwannier)
    max_at_k = [Float64[] for _ in 1:nwannier]
    total_squares = 0.0
    total_count = 0
    for iband in 1:nwannier
        errors = Float64[]
        error_indices = Int[]
        for ik in 1:nkpoints
            wannier_energy = wannier.eigenvalues[ik, iband]
            window[1] <= wannier_energy <= window[2] || continue
            qe_candidates = filter(
                energy -> window[1] <= energy <= window[2],
                view(qe_values, ik, :),
            )
            isempty(qe_candidates) && continue
            push!(errors, minimum(abs.(qe_candidates .- wannier_energy)))
            push!(error_indices, ik)
        end
        count_in_window[iband] = length(errors)
        isempty(errors) && continue

        max_index = argmax(errors)
        rms[iband] = sqrt(sum(abs2, errors) / length(errors))
        max_error[iband] = errors[max_index]
        max_at_k[iband] = copy(wannier.kpoints_frac[error_indices[max_index]])
        total_squares += sum(abs2, errors)
        total_count += length(errors)
    end

    return BandComparison(
        Matrix{Float64}(qe_values),
        count_in_window,
        rms,
        max_error,
        max_at_k,
        total_count == 0 ? NaN : sqrt(total_squares / total_count),
    )
end

function combined_rms(comparisons::BandComparison...)
    total_squares = 0.0
    total_count = 0
    for comparison in comparisons, iband in eachindex(comparison.n_in_window)
        count = comparison.n_in_window[iband]
        count == 0 && continue
        total_squares += comparison.rms_eV[iband]^2 * count
        total_count += count
    end
    return total_count == 0 ? NaN : sqrt(total_squares / total_count)
end

function write_comparison(
    path::AbstractString,
    comparison::BandComparison,
    window::Tuple{Float64, Float64},
    nkpoints::Integer,
)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# window = [$(window[1]), $(window[2])]")
        println(io, "# n_kpoints = $nkpoints")
        println(io, "# overall_rms_eV = $(comparison.overall_rms_eV)")
        println(io, "# matching = nearest QE eigenvalue at the same k-point within window")
        println(io, "# columns = wannier_band  n_in_window  rms_eV  max_eV  max_at_k")
        for iband in eachindex(comparison.n_in_window)
            kpoint = comparison.max_at_k[iband]
            ktext = isempty(kpoint) ? "nan nan nan" : join(kpoint, " ")
            @printf(
                io,
                "%d  %d  %.10f  %.10f  %s\n",
                iband,
                comparison.n_in_window[iband],
                comparison.rms_eV[iband],
                comparison.max_eV[iband],
                ktext,
            )
        end
    end
    return nothing
end

function _write_comparison_rows(io::IO, spin::Symbol, comparison::BandComparison)
    for iband in eachindex(comparison.n_in_window)
        kpoint = comparison.max_at_k[iband]
        ktext = isempty(kpoint) ? "nan nan nan" : join(kpoint, " ")
        @printf(
            io,
            "%s  %d  %d  %.10f  %.10f  %s\n",
            String(spin),
            iband,
            comparison.n_in_window[iband],
            comparison.rms_eV[iband],
            comparison.max_eV[iband],
            ktext,
        )
    end
end

function write_comparison(
    path::AbstractString,
    comparison_up::BandComparison,
    comparison_down::BandComparison,
    window::Tuple{Float64, Float64},
    nkpoints::Integer,
)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# window = [$(window[1]), $(window[2])]")
        println(io, "# n_kpoints = $nkpoints")
        println(io, "# rms_up_eV = $(comparison_up.overall_rms_eV)")
        println(io, "# rms_down_eV = $(comparison_down.overall_rms_eV)")
        println(io, "# overall_rms_eV = $(combined_rms(comparison_up, comparison_down))")
        println(io, "# matching = nearest QE eigenvalue at the same k-point within window")
        println(io, "# columns = spin  wannier_band  n_in_window  rms_eV  max_eV  max_at_k")
        _write_comparison_rows(io, :up, comparison_up)
        _write_comparison_rows(io, :down, comparison_down)
    end
    return nothing
end

end
