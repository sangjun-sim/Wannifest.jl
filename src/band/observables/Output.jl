module ObservablesOutput

using Dates
using Printf

using ..Model: BandResult, RunConfig

export oam_data_path, sam_data_path, write_oam_data, write_sam_data

function oam_data_path(bands_data_path::AbstractString)
    root, ext = splitext(String(bands_data_path))
    return string(root, "_oam", isempty(ext) ? ".dat" : ext)
end

function sam_data_path(bands_data_path::AbstractString)
    root, ext = splitext(String(bands_data_path))
    return string(root, "_sam", isempty(ext) ? ".dat" : ext)
end

function _selection_labels(config::RunConfig)
    return [string(selection.site, ":", selection.orbital_shell) for selection in config.oam.selections]
end

function write_oam_data(path::AbstractString, result::BandResult, config::RunConfig)
    oam = result.oam
    isnothing(oam) && return nothing
    mkpath(dirname(path))
    open(path, "w") do io
        nk, nbands, ncols = size(oam)
        ncols == 5 || error("OAM data has $ncols columns, expected 5")
        println(io, "# generated_at = ", string(Dates.now(Dates.UTC), "Z"))
        println(io, "# observable = projected_atomic_oam")
        println(io, "# units = hbar")
        println(io, "# l2_mode = operator_components")
        println(io, "# orbitals = ", join(_selection_labels(config), ", "))
        println(io, "# is_physical_distance = ", result.is_physical_distance)
        println(io, "# tick_positions = ", join(result.tick_positions, ", "))
        println(io, "# tick_labels = ", join(result.tick_labels, " | "))
        println(io, "# degeneracy_tol = ", config.oam.degeneracy_tol)
        println(io, "# WARNING: L2 is computed as <Lx^2 + Ly^2 + Lz^2> from PROJECTED operators.")
        println(io, "# L_norm = sqrt(<Lx>^2 + <Ly>^2 + <Lz>^2), not sqrt(<L2>).")
        println(io, "# columns = segment  kx  ky  kz  k_distance  band  energy  Lx  Ly  Lz  L_norm  L2")
        for (iseg, seg_range) in enumerate(result.segment_ranges)
            iseg > 1 && println(io)
            for ik in seg_range
                1 <= ik <= nk || error("OAM segment index $ik is out of bounds for $nk k-points")
                k = result.kpoints_frac[ik]
                for ib in 1:nbands
                    @printf(
                        io,
                        "%d  %.10f  %.10f  %.10f  %.10f  %d  %.10f  %.10f  %.10f  %.10f  %.10f  %.10f\n",
                        iseg,
                        k[1],
                        k[2],
                        k[3],
                        result.distances[ik],
                        ib,
                        result.eigenvalues[ik, ib],
                        oam[ik, ib, 1],
                        oam[ik, ib, 2],
                        oam[ik, ib, 3],
                        oam[ik, ib, 4],
                        oam[ik, ib, 5],
                    )
                end
            end
        end
    end
    return nothing
end

function write_sam_data(path::AbstractString, result::BandResult, config::RunConfig)
    sam = result.sam
    isnothing(sam) && return nothing
    mkpath(dirname(path))
    open(path, "w") do io
        nk, nbands, ncols = size(sam)
        ncols == 5 || error("SAM data has $ncols columns, expected 5")
        println(io, "# generated_at = ", string(Dates.now(Dates.UTC), "Z"))
        println(io, "# observable = spin_angular_momentum")
        println(io, "# units = hbar for Sx/Sy/Sz/S_norm, hbar^2 for S2")
        println(io, "# spin_layout = ", config.spin.layout)
        println(io, "# is_physical_distance = ", result.is_physical_distance)
        println(io, "# tick_positions = ", join(result.tick_positions, ", "))
        println(io, "# tick_labels = ", join(result.tick_labels, " | "))
        println(io, "# degeneracy_tol = ", config.sam.degeneracy_tol)
        println(io, "# S_norm = sqrt(<Sx>^2 + <Sy>^2 + <Sz>^2), not sqrt(<S2>).")
        println(io, "# columns = segment  kx  ky  kz  k_distance  band  energy  Sx  Sy  Sz  S_norm  S2")
        for (iseg, seg_range) in enumerate(result.segment_ranges)
            iseg > 1 && println(io)
            for ik in seg_range
                1 <= ik <= nk || error("SAM segment index $ik is out of bounds for $nk k-points")
                k = result.kpoints_frac[ik]
                for ib in 1:nbands
                    @printf(
                        io,
                        "%d  %.10f  %.10f  %.10f  %.10f  %d  %.10f  %.10f  %.10f  %.10f  %.10f  %.10f\n",
                        iseg,
                        k[1],
                        k[2],
                        k[3],
                        result.distances[ik],
                        ib,
                        result.eigenvalues[ik, ib],
                        sam[ik, ib, 1],
                        sam[ik, ib, 2],
                        sam[ik, ib, 3],
                        sam[ik, ib, 4],
                        sam[ik, ib, 5],
                    )
                end
            end
        end
    end
    return nothing
end

end
