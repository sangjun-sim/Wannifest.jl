module Output

using Dates
using Printf

using ..Model: BandResult, DosResult, RunConfig

export channel_path, write_bands_data, write_dos_data, write_projection_weights_data, write_pdos_data

function _metadata_path(path::AbstractString)
    return isempty(path) ? "" : abspath(path)
end

function channel_path(path::AbstractString, spin::Symbol)
    spin == :none && return String(path)
    suffix = spin == :up ? "up" : spin == :down ? "dn" : error("Unsupported channel spin: $spin")
    directory, file = dirname(path), basename(path)
    stem, extension = splitext(file)
    return joinpath(directory, "$(stem)_$(suffix)$(extension)")
end

function write_bands_data(
    path::AbstractString,
    result::BandResult;
    config::Union{Nothing, RunConfig}=nothing,
    spin::Symbol=:none,
)
    spin in (:none, :up, :down) || error("Unsupported channel spin: $spin")
    mkpath(dirname(path))
    open(path, "w") do io
        nk, nbands = size(result.eigenvalues)
        println(io, "# generated_at = ", string(Dates.now(Dates.UTC), "Z"))
        println(io, "# is_physical_distance = ", result.is_physical_distance)
        println(io, "# tick_positions = ", join(result.tick_positions, ", "))
        println(io, "# tick_labels = ", join(result.tick_labels, " | "))
        if !isnothing(config)
            println(io, "# energy_shift = ", config.energy.shift)
            spin == :none || println(io, "# spin = ", spin)
            hr_path = spin == :down ? config.files.hr_down_path : config.files.hr_path
            wsvec_path = spin == :down ? config.files.wsvec_down_path : config.files.wsvec_path
            isnothing(hr_path) && error("spin-down band data requires run.hr_down")
            println(io, "# hr = ", _metadata_path(hr_path))
            println(io, "# kpoints = ", _metadata_path(config.files.kpoints_path))
            println(io, "# structure = ", _metadata_path(config.files.structure_path))
            println(io, "# wsvec = ", isnothing(wsvec_path) ? "" : _metadata_path(wsvec_path))
        end
        print(io, "# columns = segment  kx  ky  kz  k_distance")
        for ib in 1:nbands
            print(io, "  band_", ib)
        end
        println(io)
        for (iseg, seg_range) in enumerate(result.segment_ranges)
            iseg > 1 && println(io)
            for ik in seg_range
                1 <= ik <= nk || error("Band segment index $ik is out of bounds for $nk k-points")
                k = result.kpoints_frac[ik]
                @printf(io, "%d  %.10f  %.10f  %.10f  %.10f", iseg, k[1], k[2], k[3], result.distances[ik])
                for ib in 1:nbands
                    @printf(io, "  %.10f", result.eigenvalues[ik, ib])
                end
                println(io)
            end
        end
    end
end

function write_dos_data(path::AbstractString, result::DosResult)
    mkpath(dirname(path))
    open(path, "w") do io
        if !isnothing(result.dos_down)
            println(io, "# energy  dos_up  dos_down")
            for i in eachindex(result.energies)
                @printf(io, "%.10f  %.10f  %.10f\n", result.energies[i], result.dos[i], result.dos_down[i])
            end
        else
            println(io, "# energy  dos")
            for i in eachindex(result.energies)
                @printf(io, "%.10f  %.10f\n", result.energies[i], result.dos[i])
            end
        end
    end
end

function write_projection_weights_data(path::AbstractString, result::BandResult)
    projection = result.projection
    isnothing(projection) && return nothing
    mkpath(dirname(path))
    open(path, "w") do io
        nk, nbands, ngroups = size(projection.weights)
        println(io, "# generated_at = ", string(Dates.now(Dates.UTC), "Z"))
        println(io, "# nk = $nk, nbands = $nbands, ngroups = $ngroups")
        println(io, "# labels = ", join(projection.labels, ", "))
        println(io, "# disjoint = ", projection.disjoint)
        println(io, "# covers_all = ", projection.covers_all)
        println(io, "# columns: ik ib label weight")
        for ik in 1:nk, ib in 1:nbands, ig in 1:ngroups
            @printf(io, "%d  %d  %s  %.10f\n", ik, ib, projection.labels[ig], projection.weights[ik, ib, ig])
        end
    end
    return nothing
end

function write_pdos_data(path::AbstractString, result::DosResult)
    projected = result.projected
    isnothing(projected) && return nothing
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# generated_at = ", string(Dates.now(Dates.UTC), "Z"))
        println(io, "# labels = ", join(projected.labels, ", "))
        println(io, "# atom_counts = ", join(projected.atom_counts, ", "))
        println(io, "# disjoint = ", projected.disjoint)
        println(io, "# covers_all = ", projected.covers_all)
        print(io, "# energy")
        for label in projected.labels
            print(io, "  ", label)
        end
        println(io)
        for ie in eachindex(result.energies)
            @printf(io, "%.10f", result.energies[ie])
            for ig in eachindex(projected.labels)
                @printf(io, "  %.10f", projected.pdos[ie, ig])
            end
            println(io)
        end
    end
    return nothing
end

end
