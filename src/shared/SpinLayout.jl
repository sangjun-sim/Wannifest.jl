module SpinLayout

using Printf

export DEFAULT_LAYOUT, ALLOWED_LAYOUTS, parse_layout, normalize_layout
export source_to_canonical_indices, layout_indices, band_indices
export band_color, band_plot_order, print_summary

const DEFAULT_LAYOUT = :qe
const ALLOWED_LAYOUTS = Set((:qe, :vasp544))

function _validate_layout(layout::Symbol, context::AbstractString)::Symbol
    layout in ALLOWED_LAYOUTS ||
        error("Unsupported $context: $layout (allowed: vasp544, qe)")
    return layout
end

normalize_layout(layout::Symbol; context::AbstractString="spin layout") =
    _validate_layout(layout, context)

function normalize_layout(raw::AbstractString; context::AbstractString="spin layout")::Symbol
    return _validate_layout(Symbol(lowercase(strip(String(raw)))), context)
end

normalize_layout(raw; context::AbstractString="spin layout") =
    error("$context must be a string or symbol")

function parse_layout(raw; context::AbstractString="spin.layout", default::Symbol=DEFAULT_LAYOUT)::Symbol
    isnothing(raw) && return normalize_layout(default; context=context)
    raw isa AbstractString || error("$context must be a string")
    return normalize_layout(raw; context=context)
end

function source_to_canonical_indices(num_wann::Integer, layout)::Vector{Int}
    nw = Int(num_wann)
    nw > 0 || error("spin layout requires positive num_wann, got $nw")
    mode = normalize_layout(layout; context="spin layout")
    if mode == :qe
        return collect(1:nw)
    elseif mode == :vasp544
        iseven(nw) || error("spin_layout=vasp544 requires an even num_wann, got $nw")
        half = nw ÷ 2
        # indices[source_index] = canonical_qe_index
        indices = Vector{Int}(undef, nw)
        for i in 1:half
            indices[i] = 2i - 1
            indices[i + half] = 2i
        end
        return indices
    end
    error("Unsupported spin layout: $mode")
end

function layout_indices(num_wann::Integer, layout::Symbol)
    nw = Int(num_wann)
    iseven(nw) || error("spin layout requires an even num_wann, got $nw")
    half = nw ÷ 2
    mode = normalize_layout(layout; context="spin layout")
    if mode == :qe
        return collect(1:2:nw), collect(2:2:nw)
    elseif mode == :vasp544
        return collect(1:half), collect((half + 1):nw)
    end
    error("Unsupported spin layout: $mode")
end

function band_indices(nbands::Integer, layout::Symbol)
    nb = Int(nbands)
    iseven(nb) || error("spin.enabled requires an even number of bands, got $nb")
    normalize_layout(layout; context="spin layout")
    return collect(1:2:nb), collect(2:2:nb)
end

function _is_up_index(ib::Integer, nbands::Integer, layout::Symbol)
    iseven(nbands) || error("spin.enabled requires an even number of bands, got $nbands")
    1 <= ib <= nbands || error("band index $ib is outside 1:$nbands")
    normalize_layout(layout; context="spin layout")
    return isodd(ib)
end

function band_color(config, ib::Integer, nbands::Integer)
    if config.spin.enabled
        return _is_up_index(ib, nbands, config.spin.layout) ?
            config.spin.colors[1] :
            config.spin.colors[2]
    end
    return config.band_plot.colors[mod1(ib, length(config.band_plot.colors))]
end

function band_plot_order(config, nbands::Integer)
    return collect(1:Int(nbands))
end

_format_energy(value::Float64) = isfinite(value) ? @sprintf("%.6f eV", value) : "n/a"

function _print_spin_layout_summary(io::IO, result)
    result.config.spin.enabled || return nothing
    println(io, "Spin layout: ", result.config.spin.layout)
    return nothing
end

function print_summary(result; make_plot::Bool=true, io::IO=stdout)
    println(io, "Band run complete.")
    println(io, "Hermiticity check: ", result.hermiticity_ok ? "ok" : "failed")
    _print_spin_layout_summary(io, result)
    if !isnothing(result.dos_result)
        dos = result.dos_result
        @printf(io, "DOS integral: %.6f / %d states\n", dos.integral, dos.num_bands)
        println(io, "DOS center of mass: ", _format_energy(dos.center_of_mass))
        projected = dos.projected
        if !isnothing(projected)
            println(io, "PDOS center of mass:")
            for (label, center) in zip(projected.labels, projected.centers_of_mass)
                println(io, "  ", label, ": ", _format_energy(center))
            end
        end
    end
    return nothing
end

end
