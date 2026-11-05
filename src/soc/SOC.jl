module SOC

using ..HrFormat: write_hr_blocks_raw
using ..Model: SocWannierConfig, SocCorrectionResult
using ..WannierHrIO: read_hr, pair_hermiticity_error
using ..WannierTypes: HrBlocks, RKey

export hr_path, require_hr_files, sp_hamiltonian, rotate_to_soc, soc_correction, write_soc_correction

const ROTATED_HR_NAME = "h_sp_rotated_hr.dat"
const DELTA_HR_NAME = "h_soc_delta_hr.dat"
const ORIGIN = (0, 0, 0)

hr_path(dir::AbstractString, seed::AbstractString) = joinpath(dir, seed * "_hr.dat")

_hr_paths(w::SocWannierConfig) = (up=hr_path(w.dir_up, w.seed_up), down=hr_path(w.dir_down, w.seed_down),
                                  soc=hr_path(w.dir_soc, w.seed_soc))

"""Stop before the wavefunctions are read when an `_hr.dat` of the three runs is missing."""
function require_hr_files(w::SocWannierConfig)
    for path in _hr_paths(w)
        isfile(path) || error("write_soc = true needs $path (set write_hr = .true. in all three Wannier90 runs)")
    end
    return nothing
end

# Raw blocks are subtracted R by R, which is only meaningful on one Wigner–Seitz R list.
function _require_same_R(a::HrBlocks, b::HrBlocks, label_a, label_b)
    a.normalization == b.normalization || error("$label_a and $label_b hr blocks differ in normalization")
    keys(a.hoppings) == keys(b.hoppings) ||
        error("$label_a and $label_b hr files list different R vectors; the three Wannier90 runs must share the lattice and mp_grid")
    all(a.ndegen[R] == b.ndegen[R] for R in keys(a.hoppings)) || error("$label_a and $label_b hr files differ in ndegen")
    return nothing
end

"""
    sp_hamiltonian(up, down) -> HrBlocks

`H_sp(R) = blockdiag(H_up(R), H_down(R))`: the rows of T, up block then down block,
in Wannier90 file order.
"""
function sp_hamiltonian(up::HrBlocks, down::HrBlocks)::HrBlocks
    _require_same_R(up, down, "up", "down")
    nw_up, nw = up.num_wann, up.num_wann + down.num_wann
    hoppings = Dict{RKey, Matrix{ComplexF64}}()
    for (R, H_up) in up.hoppings
        H = zeros(ComplexF64, nw, nw)
        H[1:nw_up, 1:nw_up] = H_up
        H[(nw_up + 1):nw, (nw_up + 1):nw] = down.hoppings[R]
        hoppings[R] = H
    end
    return HrBlocks("H_sp(R) = blockdiag(H_up, H_down)", nw, length(hoppings), hoppings, copy(up.ndegen), up.normalization)
end

"""
    rotate_to_soc(sp, T0) -> HrBlocks

`T(0)† H_sp(R) T(0)` for every R. With `|w^SOC_{R n}⟩ ≈ Σ_m |w^SP_{R m}⟩ T(0)_{mn}`
(the R = R′ approximation, see `weight0`) this is `⟨w^SOC_0|H_sp|w^SOC_R⟩`, so one
R-independent T serves all blocks. T(0) is used as computed, not unitarized.
"""
function rotate_to_soc(sp::HrBlocks, T0::AbstractMatrix{<:Number})::HrBlocks
    size(T0, 1) == sp.num_wann || error("T(0) has $(size(T0, 1)) rows but H_sp has $(sp.num_wann) Wannier functions")
    hoppings = Dict{RKey, Matrix{ComplexF64}}(R => T0' * H * T0 for (R, H) in sp.hoppings)
    return HrBlocks("T(0)^dagger H_sp(R) T(0), SOC Wannier basis", size(T0, 2), length(hoppings), hoppings,
                    copy(sp.ndegen), sp.normalization)
end

"""
    soc_correction(soc, up, down, T0) -> (rotated, delta)

`ΔH(R) = H_soc(R) − T(0)† H_sp(R) T(0)` in the SOC Wannier basis (columns of T),
together with the rotated SP Hamiltonian. Blocks stay raw: `ndegen` is carried over.
"""
function soc_correction(soc::HrBlocks, up::HrBlocks, down::HrBlocks, T0::AbstractMatrix{<:Number})
    size(T0, 2) == soc.num_wann || error("T(0) has $(size(T0, 2)) columns but H_soc has $(soc.num_wann) Wannier functions")
    rotated = rotate_to_soc(sp_hamiltonian(up, down), T0)
    _require_same_R(soc, rotated, "SOC", "SP")
    hoppings = Dict{RKey, Matrix{ComplexF64}}(R => H - rotated.hoppings[R] for (R, H) in soc.hoppings)
    delta = HrBlocks("H_soc(R) - T(0)^dagger H_sp(R) T(0), SOC Wannier basis", soc.num_wann, length(hoppings), hoppings,
                     copy(soc.ndegen), soc.normalization)
    return rotated, delta
end

_max_abs(hoppings, keep) = maximum((maximum(abs, H) for (R, H) in hoppings if keep(R)); init=0.0)

"""
    write_soc_correction(wannier, T0, output_dir) -> SocCorrectionResult

Read the three `_hr.dat` files in file order (no spin rearrangement, as for T),
then write `h_sp_rotated_hr.dat` and `h_soc_delta_hr.dat` in Wannier90 hr format.
"""
function write_soc_correction(w::SocWannierConfig, T0::AbstractMatrix{<:Number}, output_dir::AbstractString)::SocCorrectionResult
    paths = _hr_paths(w)
    up, down, soc = (read_hr(path; spin_layout=:qe) for path in paths)
    rotated, delta = soc_correction(soc, up, down, T0)
    haskey(delta.hoppings, ORIGIN) || error("$(paths.soc) holds no R = (0, 0, 0) block")
    save(blocks, name) = write_hr_blocks_raw(joinpath(output_dir, name), blocks.header, blocks.num_wann,
                                             blocks.hoppings, blocks.ndegen)
    return SocCorrectionResult(save(rotated, ROTATED_HR_NAME), save(delta, DELTA_HR_NAME),
                               _max_abs(delta.hoppings, ==(ORIGIN)), _max_abs(delta.hoppings, !=(ORIGIN)),
                               pair_hermiticity_error(delta))
end

end
