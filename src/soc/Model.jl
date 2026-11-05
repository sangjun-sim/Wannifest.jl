module Model

using ..WannierTypes: RKey

export SocQeConfig, SocWannierConfig, SocRunConfig, SocConfig
export QeWfcHeader, QeWavefunction, WfcFiles, WannierChannel, KpointIndex
export OnsiteDiagnostic, SocCorrectionResult, SocRunResult
export num_kpoints, effective_window

struct SocQeConfig
    sp_outdir::String
    sp_prefix::String
    sp_input::String
    soc_outdir::String
    soc_prefix::String
    soc_input::String
end

struct SocWannierConfig
    dir_up::String
    seed_up::String
    dir_down::String
    seed_down::String
    dir_soc::String
    seed_soc::String
end

struct SocRunConfig
    atol::Float64              # Cartesian k comparison, bohr⁻¹
    degeneracy_tol::Float64    # outer-window edge degeneracy warning, eV
    window_tol::Float64        # SP/SOC window comparison, eV
    eig_tol::Float64           # .eig vs QE XML eigenvalues, eV
    write_vmat::Bool
    write_soc::Bool            # ΔH(R) = H_soc − T(0)†H_sp T(0); needs the three _hr.dat
    tmat_R::Union{Symbol, Vector{RKey}}   # :hr or an explicit list of R
end

struct SocConfig
    input_path::String
    qe::SocQeConfig
    wannier::SocWannierConfig
    run::SocRunConfig
    output_dir::String
end

struct QeWfcHeader
    ik::Int
    kpoint::Vector{Float64}    # Cartesian, bohr⁻¹
    ispin::Int
    gamma_only::Bool
    scale_factor::Float64
end

struct QeWavefunction
    header::QeWfcHeader
    miller::Matrix{Int}                    # 3×nG
    coefficients::Array{ComplexF64, 3}     # nG×npol×nband
    reciprocal::Matrix{Float64}            # columns b1 b2 b3, bohr⁻¹
    ngw::Int
    stored_ng::Int                         # nG on disk, before any Gamma expansion
    source::String
end

struct WfcFiles
    up::Vector{String}
    down::Vector{String}
    soc::Vector{String}
end

"""One Wannier90 calculation: Ũ = U_dis·U with rows moved to the Wannier90 band index."""
struct WannierChannel
    label::String
    seed_path::String                      # <dir>/<seed>, without extension
    num_wann::Int
    num_bands::Int
    kpoints::Matrix{Float64}               # 3×nk fractional, `.win` order
    U::Array{ComplexF64, 3}                # num_bands×num_wann×nk
    kept::Vector{Int}                      # QE band index of every Wannier90 band
    eig::Union{Nothing, Matrix{Float64}}   # num_bands×nk, eV
    disentangled::Bool
    nfirstwin::Vector{Int}
    ndimwin::Vector{Int}
    win_min::Union{Nothing, Float64}       # dis_win_min as set in .win (eV); nothing = Wannier90 default
    win_max::Union{Nothing, Float64}
    froz_min::Union{Nothing, Float64}      # frozen window; both nothing when dis_froz_max is unset
    froz_max::Union{Nothing, Float64}
end

num_kpoints(channel::WannierChannel) = size(channel.kpoints, 2)

"""Outer window actually used by Wannier90: unset edges default to the full `.eig` range."""
function effective_window(channel::WannierChannel)
    edge(value, extremum) = !isnothing(value) ? value : isnothing(channel.eig) ? nothing : extremum(channel.eig)
    lower, upper = edge(channel.win_min, minimum), edge(channel.win_max, maximum)
    return isnothing(lower) || isnothing(upper) ? nothing : (lower, upper)
end

"""k correspondences, all indexed by the k order of the up channel."""
struct KpointIndex
    down::Vector{Int}       # k index in the down channel
    soc::Vector{Int}        # k index in the SOC channel
    sp_file::Vector{Int}    # number of the QE wfcup/wfcdw file
    soc_file::Vector{Int}   # number of the QE SOC wfc file
end

struct OnsiteDiagnostic
    weight0::Float64                       # ‖T(0)‖²_F / Σ_R ‖T(R)‖²_F
    unitarity0::Union{Nothing, Float64}    # ‖T(0)T(0)† − I‖₂
    perm0::Union{Nothing, Float64}         # ‖|T(0)| − P‖_max
end

"""`ΔH(R) = H_soc(R) − T(0)† H_sp(R) T(0)` in the SOC Wannier basis, raw Wannier90 blocks (eV)."""
struct SocCorrectionResult
    rotated_path::String                   # T(0)† H_sp(R) T(0)
    delta_path::String                     # ΔH(R)
    onsite_max_abs::Float64                # ‖ΔH(0)‖_max
    offsite_max_abs::Float64               # max_{R≠0} ‖ΔH(R)‖_max
    hermiticity::Float64                   # max_R ‖ΔH(R) − ΔH(−R)†‖_max after the ndegen division
end

struct SocRunResult
    config::SocConfig
    output_dir::String
    vmat_path::Union{Nothing, String}
    tmat_k_path::String
    tmat_r_path::String
    meta_path::String
    unitarity::Union{Nothing, Float64}     # max_k ‖T(k)T(k)† − I‖₂
    onsite::OnsiteDiagnostic
    eig_max_diff::Dict{String, Float64}    # channel label => max|ΔE| against the QE XML, eV
    window_warnings::Int
    correction::Union{Nothing, SocCorrectionResult}
end

end
