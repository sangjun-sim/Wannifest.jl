module WannierWsvecIO

using ..SpinLayout
using ..WannierTypes: RKey, WsvecEntry, WsvecTable

export read_wsvec

function _wsvec_index_map(num_wann::Union{Nothing, Int}, spin_layout)
    if isnothing(num_wann)
        SpinLayout.normalize_layout(spin_layout; context="wsvec spin_layout") == SpinLayout.DEFAULT_LAYOUT ||
            error("read_wsvec with spin_layout=vasp544 requires num_wann")
        return nothing
    end
    return SpinLayout.source_to_canonical_indices(num_wann, spin_layout)
end

function read_wsvec(
    path::AbstractString;
    num_wann::Union{Nothing, Int}=nothing,
    spin_layout=SpinLayout.DEFAULT_LAYOUT,
)::WsvecTable
    index_map = _wsvec_index_map(num_wann, spin_layout)
    open(path, "r") do io
        eof(io) && error("Empty wsvec file: $(abspath(path))")
        readline(io)

        table = Dict{Tuple{RKey, Int, Int}, WsvecEntry}()
        while !eof(io)
            line = strip(readline(io))
            isempty(line) && continue
            startswith(line, "#") && continue
            startswith(line, "!") && continue
            head = split(line)
            length(head) >= 5 || error("Malformed wsvec header line: $line")
            R = (parse(Int, head[1]), parse(Int, head[2]), parse(Int, head[3]))
            iw = parse(Int, head[4])
            jw = parse(Int, head[5])
            if !isnothing(index_map)
                1 <= iw <= length(index_map) || error("wsvec row index $iw is outside 1:$(length(index_map)) in $(abspath(path))")
                1 <= jw <= length(index_map) || error("wsvec column index $jw is outside 1:$(length(index_map)) in $(abspath(path))")
                iw = index_map[iw]
                jw = index_map[jw]
            end
            key = (R, iw, jw)
            haskey(table, key) && error("Duplicate wsvec entry for $key in $(abspath(path))")

            eof(io) && error("Unexpected EOF while reading wsvec multiplicity for $((R, iw, jw))")
            n_shift = parse(Int, strip(readline(io)))
            n_shift >= 1 || error("Invalid wsvec multiplicity for $((R, iw, jw)): n_shift must be positive")
            shifts = Matrix{Int}(undef, 3, n_shift)
            for is in 1:n_shift
                eof(io) && error("Unexpected EOF while reading wsvec shifts for $((R, iw, jw))")
                raw = split(strip(readline(io)))
                length(raw) >= 3 || error("Malformed wsvec shift line for $((R, iw, jw))")
                shifts[:, is] = parse.(Int, raw[1:3])
            end
            table[key] = WsvecEntry(n_shift, shifts)
        end
        return WsvecTable(table)
    end
end

end
