module Service

using LinearAlgebra: inv
using TOML

using ..BandWindow: check_eig_against_qe
using ..BinaryIO: write_tensor, HEADER_BYTES
using ..Channel: load_channel, warn_window_mismatch, warn_projection_mismatch
using ..KpointMatch: match_kpoints, require_same_kpoints
using ..Model
using ..QeInputCheck: check_qe_inputs
using ..QeWfcIO: qe_save_dir, sp_wfc_files, soc_wfc_files
using ..QeXmlSchema: qe_wfc_metadata, qe_eigenvalues
using ..SOC: hr_path, require_hr_files, write_soc_correction
using ..Tmat: tmat_k, tmat_r, tmat_unitarity_error, tmat_onsite_diagnostic, warn_onsite_diagnostic
using ..Validate: validate_same_lattice
using ..Vmat: wfc_kpoints, calculate_vmat
using ..WannierHrIO
using ..WannierTypes: RKey

export run

const ROW_COLUMN_RULE = "row i (up, i = 1..nw_up) <-> SOC column 2i-1; row nw_up+i (down, i = 1..nw_down) <-> SOC column 2i; " *
    "expected T(0) ≈ P, the permutation of this rule (Wannier90 v3 spinor order), when nw_up == nw_down and nw_soc == 2 nw_up"
const DELTA_HR_CONVENTION = "ΔH(R) = H_soc(R) − T(0)† H_sp(R) T(0) with H_sp = blockdiag(H_up, H_down), SOC Wannier basis in file order; " *
    "raw Wannier90 hr blocks (eV) on the R list and ndegen of the input hr files"
const TMAT_R_CONVENTION = "T(R) = (1/N) sum_k exp(+2πi k·R) T(k) = <w^SP_R|w^SOC_0> = <w^SP_0|w^SOC_{-R}>. " *
    "Wannier90 hr uses H(R) = (1/N) sum_k exp(-2πi k·R) H(k) = <w_0|H|w_R>, so the hr-convention matrix is T_hr(R) = <w^SP_0|w^SOC_R> = T(-R)"

function _load_channels(config::SocConfig, sp_meta, soc_meta)
    w, run = config.wannier, config.run
    load(dir, seed, label, nbnd_qe) = load_channel(dir, seed; label, nbnd_qe, atol=run.atol, degeneracy_tol=run.degeneracy_tol)
    return load(w.dir_up, w.seed_up, "up", sp_meta.nbnd), load(w.dir_down, w.seed_down, "down", sp_meta.nbnd_down),
           load(w.dir_soc, w.seed_soc, "soc", soc_meta.nbnd)
end

"""k correspondences of the down/SOC channels and of the QE files, all in up-channel k order."""
function _kpoint_index(up, down, soc, reciprocal, sp_kcart, soc_kcart, atol)
    requested = reciprocal * up.kpoints
    return KpointIndex(
        match_kpoints(requested, reciprocal * down.kpoints; atol, label="down Wannier channel"),
        match_kpoints(requested, reciprocal * soc.kpoints; atol, label="SOC Wannier channel"),
        match_kpoints(requested, sp_kcart; atol, label="SP QE wavefunctions"),
        match_kpoints(requested, soc_kcart; atol, label="SOC QE wavefunctions"),
    )
end

# QE k index for the channel's own k order: own[position[k]] = file[k].
function _own_order(position::Vector{Int}, file::Vector{Int})
    own = similar(file)
    own[position] = file
    return own
end

"""§1.2.1: `.eig` of every channel against the QE XML eigenvalues, before V is computed."""
function _check_eigenvalues(up, down, soc, index, sp_save, soc_save, sp_kcart, soc_kcart, run)
    qe_sp, qe_soc = qe_eigenvalues(sp_save; soc=false), qe_eigenvalues(soc_save; soc=true)
    require_same_kpoints(sp_kcart, qe_sp.kcart; atol=run.atol, label="SP data-file-schema.xml")
    require_same_kpoints(soc_kcart, qe_soc.kcart; atol=run.atol, label="SOC data-file-schema.xml")
    cases = ((up, qe_sp.up, index.sp_file), (down, qe_sp.down, _own_order(index.down, index.sp_file)),
             (soc, qe_soc.up, _own_order(index.soc, index.soc_file)))
    max_diff = Dict{String, Float64}()
    for (channel, qe_eig, kindex) in cases
        isnothing(channel.eig) && continue
        max_diff[channel.label] = check_eig_against_qe(channel.eig, qe_eig, channel.kept, kindex;
                                                       eig_tol=run.eig_tol, label=channel.label)
    end
    return max_diff
end

function _tmat_R(config::SocConfig)
    config.run.tmat_R isa Symbol || return config.run.tmat_R, nothing
    path = hr_path(config.wannier.dir_up, config.wannier.seed_up)
    isfile(path) || error("tmat_R = \"hr\" needs $path (set write_hr = .true. in the up run)")
    hr = WannierHrIO.read_hr(path; spin_layout=:qe)
    Rs = sort!(collect(keys(hr.ndegen)))
    return Rs, [hr.ndegen[R] for R in Rs]
end

_R_labels(Rs::Vector{RKey}) = Int64[R[d] for d in 1:3, R in Rs]

function _window_summary(channel::WannierChannel)
    summary = Dict{String, Any}(
        "num_bands" => channel.num_bands, "num_wann" => channel.num_wann, "disentangled" => channel.disentangled,
        "nfirstwin_min" => minimum(channel.nfirstwin), "nfirstwin_max" => maximum(channel.nfirstwin),
        "ndimwin_min" => minimum(channel.ndimwin), "ndimwin_max" => maximum(channel.ndimwin),
        "kpoints_with_narrowed_window" => count(<(channel.num_bands), channel.ndimwin),
    )
    window = effective_window(channel)
    isnothing(window) || (summary["outer_window_eV"] = collect(window))
    isnothing(channel.froz_max) || (summary["frozen_window_eV"] = [channel.froz_min, channel.froz_max])
    return summary
end

_set_optional!(table, key, value) = isnothing(value) ? table : setindex!(table, value, key)

function _write_meta(path, config, channels, dims, outputs, unitarity, onsite, eig_max_diff, Rs, ndegen, correction)
    up, down, soc = channels
    diagnostics = Dict{String, Any}("weight0" => onsite.weight0, "eig_max_abs_diff_eV" => eig_max_diff,
                                    "unitarity_norm" => "spectral norm of T T† − I")
    _set_optional!(diagnostics, "unitarity_max_k", unitarity)
    _set_optional!(diagnostics, "unitarity0", onsite.unitarity0)
    _set_optional!(diagnostics, "perm0", onsite.perm0)
    if !isnothing(correction)
        diagnostics["delta_hr"] = Dict{String, Any}("onsite_max_abs_eV" => correction.onsite_max_abs,
                                                    "offsite_max_abs_eV" => correction.offsite_max_abs,
                                                    "hermiticity_error_eV" => correction.hermiticity)
    end
    meta = Dict{String, Any}(
        "format" => Dict{String, Any}("magic" => "WNFSOC01", "header_bytes" => HEADER_BYTES, "byte_order" => "little",
                                      "data" => "ComplexF64 [n1, n2, nblock], column-major, after 3×nblock block labels"),
        "dimensions" => dims,
        "layout" => Dict{String, Any}(
            "k_order" => "k-points of the up-channel .win (fractional), stored as block labels of vmat and tmat_k",
            "vmat_rows" => "kept SP bands: up block (kept_bands.up), then down block (kept_bands.down)",
            "vmat_columns" => "kept SOC bands (kept_bands.soc)",
            "tmat_rows" => "SP Wannier functions: up block, then down block, in Wannier90 file order",
            "tmat_columns" => "SOC Wannier functions in Wannier90 file order (columns are not rearranged)",
            "tmat_row_column_rule" => ROW_COLUMN_RULE,
        ),
        "tmat_r_convention" => TMAT_R_CONVENTION,
        "tmat_r_R" => [collect(R) for R in Rs],
        "kept_bands" => Dict{String, Any}(channel.label => channel.kept for channel in channels),
        "window" => Dict{String, Any}(channel.label => _window_summary(channel) for channel in channels),
        "run" => Dict{String, Any}("atol" => config.run.atol, "degeneracy_tol" => config.run.degeneracy_tol,
                                   "window_tol" => config.run.window_tol, "eig_tol" => config.run.eig_tol,
                                   "tmat_R" => config.run.tmat_R isa Symbol ? "hr" : "list"),
        "inputs" => Dict{String, Any}(
            "input" => config.input_path, "sp_save" => qe_save_dir(config.qe.sp_outdir, config.qe.sp_prefix),
            "soc_save" => qe_save_dir(config.qe.soc_outdir, config.qe.soc_prefix),
            "sp_input" => config.qe.sp_input, "soc_input" => config.qe.soc_input,
            "seed_up" => up.seed_path, "seed_down" => down.seed_path, "seed_soc" => soc.seed_path,
        ),
        "outputs" => outputs,
        "diagnostics" => diagnostics,
    )
    # Raw T(R) repeats on Wigner–Seitz boundary R; divide by ndegen(R) before summing over R.
    _set_optional!(meta, "tmat_r_ndegen", ndegen)
    _set_optional!(meta, "delta_hr_convention", isnothing(correction) ? nothing : DELTA_HR_CONVENTION)
    open(io -> TOML.print(io, meta; sorted=true), path, "w")
    return path
end

"""
    run(config; write_vmat=config.run.write_vmat, output_dir=config.output_dir, progress=false) -> SocRunResult

Check the QE inputs, load the three Wannier channels, match all k lists, verify
`.eig` against the QE XML, then stream V(k) and write `vmat.bin` (optional),
`tmat_k.bin`, `tmat_r.bin` and `soc_meta.toml` into `output_dir`. With
`write_soc`, `T(0)` also gives `h_sp_rotated_hr.dat` and `h_soc_delta_hr.dat`.
"""
function run(config::SocConfig; write_vmat::Bool=config.run.write_vmat,
             output_dir::AbstractString=config.output_dir, progress::Bool=false)::SocRunResult
    config.run.write_soc && require_hr_files(config.wannier)
    check_qe_inputs(config.qe.sp_input, config.qe.soc_input)
    sp_save = qe_save_dir(config.qe.sp_outdir, config.qe.sp_prefix)
    soc_save = qe_save_dir(config.qe.soc_outdir, config.qe.soc_prefix)
    sp_meta, soc_meta = qe_wfc_metadata(sp_save; soc=false), qe_wfc_metadata(soc_save; soc=true)
    validate_same_lattice(sp_meta.cell, soc_meta.cell; atol=config.run.atol)

    channels = _load_channels(config, sp_meta, soc_meta)
    up, down, soc = channels
    window_warnings = warn_window_mismatch(up, down, soc; window_tol=config.run.window_tol) +
        warn_projection_mismatch((channel.seed_path * ".win" for channel in channels)...)

    files = WfcFiles(sp_wfc_files(sp_save, sp_meta.nk)..., soc_wfc_files(soc_save, soc_meta.nk))
    sp_kcart, soc_kcart = wfc_kpoints(files.up), wfc_kpoints(files.soc)
    reciprocal = 2pi .* inv(sp_meta.cell)'
    index = _kpoint_index(up, down, soc, reciprocal, sp_kcart, soc_kcart, config.run.atol)
    eig_max_diff = _check_eigenvalues(up, down, soc, index, sp_save, soc_save, sp_kcart, soc_kcart, config.run)

    V = calculate_vmat(files, sp_meta, soc_meta, index, up.kept, down.kept, soc.kept; progress)
    mkpath(output_dir)
    vmat_path = write_vmat ? write_tensor(joinpath(output_dir, "vmat.bin"), "vmat", V, up.kpoints) : nothing
    Tk = tmat_k(up, down, soc, V, index)
    tmat_k_path = write_tensor(joinpath(output_dir, "tmat_k.bin"), "tmat_k", Tk, up.kpoints)
    Rs, ndegen = _tmat_R(config)
    tmat_r_path = write_tensor(joinpath(output_dir, "tmat_r.bin"), "tmat_r", tmat_r(Tk, up.kpoints, Rs), _R_labels(Rs))

    T0 = tmat_r(Tk, up.kpoints, RKey[(0, 0, 0)])[:, :, 1]
    onsite = tmat_onsite_diagnostic(Tk, T0, up.num_wann, down.num_wann)
    warn_onsite_diagnostic(onsite)
    unitarity = tmat_unitarity_error(Tk)
    correction = config.run.write_soc ? write_soc_correction(config.wannier, T0, output_dir) : nothing

    dims = Dict{String, Any}("nk" => size(V, 3), "nR" => length(Rs), "nb_up" => up.num_bands, "nb_down" => down.num_bands,
                             "nb_soc" => soc.num_bands, "nw_up" => up.num_wann, "nw_down" => down.num_wann, "nw_soc" => soc.num_wann)
    outputs = Dict{String, Any}("tmat_k" => basename(tmat_k_path), "tmat_r" => basename(tmat_r_path))
    _set_optional!(outputs, "vmat", isnothing(vmat_path) ? nothing : basename(vmat_path))
    if !isnothing(correction)
        outputs["h_sp_rotated_hr"], outputs["h_soc_delta_hr"] = basename(correction.rotated_path), basename(correction.delta_path)
    end
    meta_path = _write_meta(joinpath(output_dir, "soc_meta.toml"), config, channels, dims, outputs,
                            unitarity, onsite, eig_max_diff, Rs, ndegen, correction)
    return SocRunResult(config, String(output_dir), vmat_path, tmat_k_path, tmat_r_path, meta_path,
                        unitarity, onsite, eig_max_diff, window_warnings, correction)
end

end
