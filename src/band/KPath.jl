module KPath

using LinearAlgebra

using ..LatticeIO: LatticeData

export KSegment, KPathData, KPathResult, parse_kpoints, generate_kpath

struct KSegment
    k_start::Vector{Float64}
    k_end::Vector{Float64}
    label_start::String
    label_end::String
end

struct KPathData
    segments::Vector{KSegment}
    npts_per_segment::Int
end

struct KPathResult
    kpoints::Vector{Vector{Float64}}
    distances::Vector{Float64}
    tick_positions::Vector{Float64}
    tick_labels::Vector{String}
    segment_ranges::Vector{UnitRange{Int}}
    is_physical_distance::Bool
end

function Base.iterate(result::KPathResult, state::Int=1)
    state == 1 && return (result.kpoints, 2)
    state == 2 && return (result.distances, 3)
    state == 3 && return (result.tick_positions, 4)
    state == 4 && return (result.tick_labels, 5)
    state == 5 && return (result.segment_ranges, 6)
    state == 6 && return (result.is_physical_distance, 7)
    return nothing
end

Base.length(::KPathResult) = 6

function parse_kpoints(path::AbstractString; lattice::Union{Nothing, LatticeData}=nothing)
    lines = readlines(path)
    length(lines) >= 4 || error("KPOINTS file too short: $path")

    npts = parse(Int, strip(lines[2]))
    npts >= 2 || error("KPOINTS line-mode requires at least 2 points per segment")
    mode_line = lowercase(strip(lines[3]))
    startswith(mode_line, "l") || error("Only line-mode KPOINTS supported, got: $(lines[3])")
    coord_line = lowercase(strip(lines[4]))

    is_cartesian = startswith(coord_line, "c")
    is_reciprocal = startswith(coord_line, "r")
    (!is_cartesian && !is_reciprocal) && error("KPOINTS header must be reciprocal or cartesian, got: $(lines[4])")
    if is_cartesian && isnothing(lattice)
        error("Cartesian KPOINTS requires a lattice file. Provide [run] structure as POSCAR/CONTCAR or wannier90.win.")
    end

    segments = KSegment[]
    idx = 5
    while idx <= length(lines)
        while idx <= length(lines) && isempty(strip(lines[idx]))
            idx += 1
        end
        idx > length(lines) && break
        first_line = idx
        k1, label1 = _parse_kpoint_line(lines[idx], path, idx)
        idx += 1

        while idx <= length(lines) && isempty(strip(lines[idx]))
            idx += 1
        end
        idx > length(lines) && error("Incomplete k-point pair in KPOINTS after line $first_line in $path")
        k2, label2 = _parse_kpoint_line(lines[idx], path, idx)
        idx += 1

        if is_cartesian
            k1 = _cart_to_reduced(k1, lattice)
            k2 = _cart_to_reduced(k2, lattice)
        end
        _same_kpoint(k1, k2) && error("Zero-length KPOINTS segment starting at line $first_line in $path")
        push!(segments, KSegment(k1, k2, label1, label2))
    end

    isempty(segments) && error("No line-mode segments found in KPOINTS: $path")
    return KPathData(segments, npts)
end

function _parse_kpoint_line(line::AbstractString, path::AbstractString, line_number::Integer)
    clean = strip(line)
    parts = split(clean)
    length(parts) >= 3 || error("Malformed KPOINTS line $line_number in $path: $line")
    k = try
        [parse(Float64, parts[1]), parse(Float64, parts[2]), parse(Float64, parts[3])]
    catch
        error("Malformed KPOINTS numeric value on line $line_number in $path: $line")
    end
    bang = findfirst('!', line)
    label = isnothing(bang) ? "" : strip(line[nextind(line, bang):end])
    return k, label
end

function _cart_to_reduced(k_cart::Vector{Float64}, lattice::LatticeData)
    return vec(lattice.reciprocal_lattice \ k_cart)
end

function generate_kpath(kpd::KPathData; lattice::Union{Nothing, LatticeData}=nothing)
    kpoints = Vector{Float64}[]
    distances = Float64[]
    tick_pos = Float64[]
    tick_labels = String[]
    segment_ranges = UnitRange{Int}[]
    current_dist = 0.0

    for (iseg, seg) in enumerate(kpd.segments)
        npts = kpd.npts_per_segment
        if iseg == 1
            push!(kpoints, copy(seg.k_start))
            push!(distances, current_dist)
            push!(tick_pos, current_dist)
            push!(tick_labels, _format_tick_label(seg.label_start; reduced_mode=isnothing(lattice)))
        elseif _same_kpoint(kpd.segments[iseg - 1].k_end, seg.k_start)
            tick_labels[end] = _merge_tick_labels(
                kpd.segments[iseg - 1].label_end,
                seg.label_start;
                reduced_mode=isnothing(lattice),
            )
        else
            tick_labels[end] = _merge_tick_labels(
                kpd.segments[iseg - 1].label_end,
                seg.label_start;
                reduced_mode=isnothing(lattice),
            )
            push!(kpoints, copy(seg.k_start))
            push!(distances, current_dist)
        end

        segment_start = length(kpoints)
        for i in 1:(npts - 1)
            t = i / (npts - 1)
            k = (1 - t) .* seg.k_start .+ t .* seg.k_end
            current_dist += _segment_step(k .- kpoints[end], lattice)
            push!(kpoints, k)
            push!(distances, current_dist)
        end

        push!(segment_ranges, segment_start:length(kpoints))
        push!(tick_pos, current_dist)
        push!(tick_labels, _format_tick_label(seg.label_end; reduced_mode=isnothing(lattice)))
    end

    return KPathResult(kpoints, distances, tick_pos, tick_labels, segment_ranges, !isnothing(lattice))
end

function _segment_step(dk::Vector{Float64}, lattice::Union{Nothing, LatticeData})
    if isnothing(lattice)
        return norm(dk)
    end
    return norm(lattice.reciprocal_lattice * dk)
end

function _same_kpoint(k1::Vector{Float64}, k2::Vector{Float64})
    return all(isapprox.(k1, k2; atol=1e-10, rtol=0.0))
end

function _merge_tick_labels(left::String, right::String; reduced_mode::Bool)
    parts = filter(!isempty, [_format_point_label(left), _format_point_label(right)])
    unique!(parts)
    merged = join(parts, "|")
    if reduced_mode && !isempty(merged)
        return string(merged, " (reduced)")
    end
    return merged
end

function _format_tick_label(label::String; reduced_mode::Bool)
    text = _format_point_label(label)
    if reduced_mode && !isempty(text)
        return string(text, " (reduced)")
    end
    return text
end

function _format_point_label(label::String)
    text = strip(label)
    isempty(text) && return ""
    text = replace(text, raw"\Gamma" => "Γ")

    io = IOBuffer()
    i = firstindex(text)
    while i <= lastindex(text)
        if text[i] == '_'
            j = nextind(text, i)
            if j <= lastindex(text) && isdigit(text[j])
                while j <= lastindex(text) && isdigit(text[j])
                    print(io, _subscript_char(text[j]))
                    j = nextind(text, j)
                end
                i = j
                continue
            end
        end
        print(io, text[i])
        i = nextind(text, i)
    end
    return String(take!(io))
end

function _subscript_char(c::Char)
    return if c == '0'
        '₀'
    elseif c == '1'
        '₁'
    elseif c == '2'
        '₂'
    elseif c == '3'
        '₃'
    elseif c == '4'
        '₄'
    elseif c == '5'
        '₅'
    elseif c == '6'
        '₆'
    elseif c == '7'
        '₇'
    elseif c == '8'
        '₈'
    elseif c == '9'
        '₉'
    else
        c
    end
end

end
