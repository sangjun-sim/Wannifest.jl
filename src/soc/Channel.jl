module Channel

using ..BandWindow
using ..KpointMatch: require_same_kpoints
using ..Model: WannierChannel, effective_window
using ..WannierEigIO: read_eig
using ..WannierUmatIO: read_umat_file
using ..Win90Keywords
using ..Win90ProjectionIO: read_win_projection_source

export load_channel, warn_window_mismatch, warn_projection_mismatch

const WIN_KEYS = Set(["num_wann", "num_bands", "write_u_matrices", "exclude_bands", "dis_spheres_num",
                      "dis_win_min", "dis_win_max", "dis_froz_min", "dis_froz_max"])
const OFFDIAGONAL_TOL = 1e-6    # eV, U_dis† diag(E) U_dis must be diagonal (text precision is 1e-10)
const CENTER_TOL = 1e-3         # fractional, projection site comparison

function _win_kpoints(win::AbstractString)
    lines = read_win_block(win, "kpoints")
    (isnothing(lines) || isempty(lines)) && error("Missing or empty kpoints block in $win")
    kpoints = Matrix{Float64}(undef, 3, length(lines))
    for (k, line) in enumerate(lines)
        fields = split(line)
        length(fields) >= 3 || error("$win: kpoints line $k needs three coordinates")
        kpoints[:, k] = parse_win_real.(fields[1:3])
    end
    return kpoints
end

_optional_real(keys, name) = haskey(keys, name) ? parse_win_real(keys[name]) : nothing

function _dimensions(keys, win)
    haskey(keys, "write_u_matrices") && parse_win_bool(keys["write_u_matrices"]) ||
        error("$win: set write_u_matrices = .true. so that Wannier90 writes <seed>_u.mat and <seed>_u_dis.mat")
    haskey(keys, "num_wann") || error("$win: missing num_wann")
    num_wann = parse_win_int(keys["num_wann"])
    num_bands = haskey(keys, "num_bands") ? parse_win_int(keys["num_bands"]) : num_wann
    num_wann > 0 || error("$win: invalid num_wann = $num_wann")
    num_bands >= num_wann || error("$win: num_bands = $num_bands is smaller than num_wann = $num_wann")
    parse_win_int(get(keys, "dis_spheres_num", "0")) > 0 &&
        error("$win: dis_spheres_num > 0 changes the outer-window rule and is not supported")
    return num_wann, num_bands
end

function _warn_window_quality(label, u_dis, eig, nfirstwin, ndimwin, degeneracy_tol)
    offdiagonal, k = subspace_offdiagonal(u_dis, eig, nfirstwin, ndimwin)
    offdiagonal <= OFFDIAGONAL_TOL ||
        @warn "$label: U_dis† diag(E) U_dis is not diagonal (max off-diagonal $offdiagonal eV at k-point $k). The rebuilt outer window may be shifted against the Wannier90 run (stale .eig or dis_win_*), unless site symmetry was used."
    flagged = window_edge_degeneracies(eig, nfirstwin, ndimwin; degeneracy_tol)
    isempty(flagged) ||
        @warn "$label: the outer window cuts a (near-)degenerate multiplet at $(length(flagged)) k-point(s), first at k-point $(first(flagged)) (|ΔE| < $degeneracy_tol eV); the subspace is ambiguous there."
    return nothing
end

"""
    load_channel(dir, seed; label, nbnd_qe, atol, degeneracy_tol) -> WannierChannel

Read `<seed>.win`, `<seed>.eig`, `<seed>_u.mat` and, for a disentangled run,
`<seed>_u_dis.mat`, and build Ũ = U_dis·U on Wannier90 band indices. Without
disentanglement (`num_bands == num_wann`) Ũ = U. The `.win` k list and the k lists
of both matrix files must agree in order, because the outer window is rebuilt
per k index from `.eig`.
"""
function load_channel(dir::AbstractString, seed::AbstractString; label::AbstractString, nbnd_qe::Integer,
                      atol::Real, degeneracy_tol::Real)::WannierChannel
    seed_path = joinpath(dir, seed)
    win = seed_path * ".win"
    keys = read_win_keywords(win, WIN_KEYS)
    num_wann, num_bands = _dimensions(keys, win)
    exclude = haskey(keys, "exclude_bands") ? parse_band_ranges(keys["exclude_bands"]) : Int[]
    kept = kept_bands(nbnd_qe, exclude, num_bands)
    kpoints = _win_kpoints(win)
    nk = size(kpoints, 2)
    u = read_umat_file(seed_path * "_u.mat"; expected_rows=num_wann, expected_cols=num_wann)
    require_same_kpoints(kpoints, u.kpoints; atol, label="$(seed)_u.mat")

    disentangled = num_bands > num_wann
    win_min = disentangled ? _optional_real(keys, "dis_win_min") : nothing
    win_max = disentangled ? _optional_real(keys, "dis_win_max") : nothing
    froz_max = disentangled ? _optional_real(keys, "dis_froz_max") : nothing
    froz_min = isnothing(froz_max) ? nothing : something(_optional_real(keys, "dis_froz_min"), -Inf)
    eig_path = seed_path * ".eig"
    eig = isfile(eig_path) ? read_eig(eig_path, num_bands, nk) : nothing
    if isnothing(eig)
        isnothing(win_min) && isnothing(win_max) ||
            error("$eig_path not found, but dis_win_min/max is set: the outer window cannot be rebuilt")
        @warn "$label: $eig_path not found; the .eig ↔ QE XML eigenvalue cross-check is skipped."
    end

    nfirstwin, ndimwin = ones(Int, nk), fill(num_bands, nk)
    if disentangled
        u_dis = read_umat_file(seed_path * "_u_dis.mat"; expected_rows=num_bands, expected_cols=num_wann)
        require_same_kpoints(kpoints, u_dis.kpoints; atol, label="$(seed)_u_dis.mat")
        if !isnothing(eig)
            nfirstwin, ndimwin = outer_window(eig, something(win_min, minimum(eig)), something(win_max, maximum(eig)), num_wann)
        end
        U = aligned_u(u_dis, u, nfirstwin, ndimwin)
        isnothing(eig) || _warn_window_quality(label, u_dis, eig, nfirstwin, ndimwin, degeneracy_tol)
    else
        isfile(seed_path * "_u_dis.mat") &&
            @warn "$label: num_bands == num_wann, so $(seed)_u_dis.mat is ignored as a leftover of another run."
        U = copy(u.data)
    end
    return WannierChannel(label, seed_path, num_wann, num_bands, kpoints, U, kept, eig, disentangled,
                          nfirstwin, ndimwin, win_min, win_max, froz_min, froz_max)
end

_differs(a, b, tol) = a != b && !(abs(a - b) <= tol)    # equal infinities agree

function _sp_union(up, down)
    isnothing(up) && return down
    isnothing(down) && return up
    return (min(up[1], down[1]), max(up[2], down[2]))
end

function _warn_window(kind, keyword, soc, expected, tol)
    isnothing(soc) && isnothing(expected) && return 0
    if isnothing(soc) || isnothing(expected) || _differs(soc[1], expected[1], tol) || _differs(soc[2], expected[2], tol)
        describe(window) = isnothing(window) ? "(none)" : "[$(window[1]), $(window[2])] eV"
        @warn "SOC $kind window $(describe(soc)) differs from the union of SP up/down windows $(describe(expected)); set $keyword of the SOC run to the union so that T(k) stays close to unitary."
        return 1
    end
    return 0
end

"""
    warn_window_mismatch(up, down, soc; window_tol) -> Int

T is close to unitary only when the SP and SOC subspaces cover the same energy
range. Warn (never stop) when the SOC outer or frozen window differs from the
union of the SP windows, or when `nw_soc ≠ nw_up + nw_down`; return the number of
warnings. Unset edges count as the full `.eig` range; an edge left at its default
in all three runs is skipped, since every channel keeps all of its bands there.
"""
function warn_window_mismatch(up::WannierChannel, down::WannierChannel, soc::WannierChannel; window_tol::Real)::Int
    count = 0
    if soc.num_wann != up.num_wann + down.num_wann
        @warn "nw_soc = $(soc.num_wann) differs from nw_up + nw_down = $(up.num_wann + down.num_wann); T(k) is rectangular and cannot be unitary."
        count += 1
    end
    windows = map(effective_window, (up, down, soc))
    if !any(isnothing, windows)
        expected, actual = collect(_sp_union(windows[1], windows[2])), collect(windows[3])
        for (edge, field) in enumerate((:win_min, :win_max))
            # All defaults: nothing to align on this side.
            all(channel -> isnothing(getfield(channel, field)), (up, down, soc)) && (actual[edge] = expected[edge])
        end
        count += _warn_window("outer", "dis_win_min/max", Tuple(actual), Tuple(expected), window_tol)
    end
    frozen(channel) = isnothing(channel.froz_max) ? nothing : (channel.froz_min, channel.froz_max)
    count += _warn_window("frozen", "dis_froz_min/max", frozen(soc), _sp_union(frozen(up), frozen(down)), window_tol)
    return count
end

function _projection_signature(win::AbstractString)
    source = read_win_projection_source(win)
    seeds = [(seed.site_label, seed.orbital, seed.center_frac) for seed in source.seeds]
    # Text after the second colon holds the local axes (z=..., x=...) of each projection line.
    axes = [lowercase(replace(get(split(line, ':'; limit=3), 3, ""), r"\s+" => ""))
            for line in source.metadata["raw_projection_lines"]]
    return seeds, axes
end

function _same_projections(a, b)
    (seeds_a, axes_a), (seeds_b, axes_b) = a, b
    length(seeds_a) == length(seeds_b) && axes_a == axes_b || return false
    return all(zip(seeds_a, seeds_b)) do (x, y)
        shift = collect(x[3]) .- collect(y[3])
        x[1] == y[1] && x[2] == y[2] && maximum(abs, shift .- round.(shift)) <= CENTER_TOL
    end
end

"""
    warn_projection_mismatch(win_up, win_down, win_soc) -> Int

T(0) is close to the permutation P only when the three `.win` files list the same
projections (sites, orbital order, local axes). Warn when they differ or cannot
be compared; return the number of warnings.
"""
function warn_projection_mismatch(win_up::AbstractString, win_down::AbstractString, win_soc::AbstractString)::Int
    signatures = try
        map(_projection_signature, (win_up, win_down, win_soc))
    catch err
        @warn "Could not compare the projections of the three .win files: $(sprint(showerror, err))"
        return 1
    end
    count = 0
    for (name, signature) in (("down", signatures[2]), ("SOC", signatures[3]))
        _same_projections(signatures[1], signature) && continue
        @warn "The $name .win lists projections that differ from the up .win (site, orbital order or local axes); T(0) will not be close to the permutation P."
        count += 1
    end
    return count
end

end
