module Win90Basis

using ..OrbitalProjection: ProjectionGroup, ProjectionSpec
using ..SpinLayout
using ..Win90OrbitalTokens: expand_orbital_token, expand_orbitals_and_shells
using ..Win90ProjectionIO: WinProjectionSeed, read_win_projection_source

export WinOrbital, WinBasis, read_win_basis, build_projection_spec
export expand_orbital_token, expand_orbitals_and_shells

struct WinOrbital
    index::Int
    species::String
    site::Int
    site_label::String
    orbital::String
    spin::Symbol
    center_frac::NTuple{3, Float64}
end

struct WinBasis
    num_wann::Int
    spinors::Bool
    species_atoms::Dict{String, Vector{NTuple{3, Float64}}}
    orbitals::Vector{WinOrbital}
    metadata::Dict{String, Any}
end

function _push_win_orbital!(out::Vector{WinOrbital}, seed::WinProjectionSeed, spin::Symbol)
    push!(out, WinOrbital(length(out) + 1, seed.species, seed.site, seed.site_label, seed.orbital, spin, seed.center_frac))
    return out
end

function _materialize_orbitals(
    seeds::Vector{WinProjectionSeed},
    spinors::Bool,
    layout::Symbol,
)
    orbitals = WinOrbital[]
    if !spinors
        for seed in seeds
            _push_win_orbital!(orbitals, seed, :unpolarized)
        end
        return orbitals
    end

    if layout == :qe
        for seed in seeds
            _push_win_orbital!(orbitals, seed, :up)
            _push_win_orbital!(orbitals, seed, :dn)
        end
    elseif layout == :vasp544
        for spin in (:up, :dn)
            for seed in seeds
                _push_win_orbital!(orbitals, seed, spin)
            end
        end
    else
        error("Unsupported spin layout: $layout")
    end
    return orbitals
end

function read_win_basis(path::AbstractString; spin_layout=SpinLayout.DEFAULT_LAYOUT)
    layout_mode = SpinLayout.normalize_layout(spin_layout; context="win spin_layout")
    source = read_win_projection_source(path)
    metadata = copy(source.metadata)
    metadata["spin_layout"] = String(layout_mode)

    orbitals = _materialize_orbitals(source.seeds, source.spinors, layout_mode)
    length(orbitals) == source.num_wann ||
        error("win basis says num_wann=$(source.num_wann), projections expand to $(length(orbitals))")
    seen = Set{Tuple{String, String, Symbol}}()
    for orbital in orbitals
        key = (orbital.site_label, orbital.orbital, orbital.spin)
        key in seen && error("duplicate win orbital entry $(key)")
        push!(seen, key)
    end
    if source.spinors
        count(o -> o.spin == :up, orbitals) == count(o -> o.spin == :dn, orbitals) ||
            error("spinors=.true. requires matching up/dn orbital counts")
    end

    return WinBasis(source.num_wann, source.spinors, source.species_atoms, orbitals, metadata)
end

function _group_has_orbital_selector(group)
    return !isempty(group.orbitals) || !isempty(group.orbital_shells)
end

function _group_has_site_selector(group)
    return !isempty(group.species) || !isempty(group.sites) || !isempty(group.site_labels)
end

function _matches_site(group, orbital::WinOrbital)
    species_ok = isempty(group.species) || orbital.species in group.species
    site_ok = if isempty(group.sites) && isempty(group.site_labels)
        true
    else
        (orbital.site in group.sites) ||
            (orbital.site_label in group.site_labels)
    end
    return species_ok && site_ok
end

function _matches_spin(group, orbital::WinOrbital)
    spin = lowercase(strip(group.spin))
    (isempty(spin) || spin == "any") && return true
    spin == "up" && return orbital.spin == :up
    spin == "dn" && return orbital.spin == :dn
    spin == "down" && return orbital.spin == :dn
    error("group '$(group.label)': spin must be 'up', 'dn', or 'any'")
end

function _selector_description(group)
    parts = String[]
    !isempty(group.species) && push!(parts, "species=$(join(group.species, ","))")
    !isempty(group.sites) && push!(parts, "sites=$(join(group.sites, ","))")
    !isempty(group.site_labels) && push!(parts, "site_labels=$(join(group.site_labels, ","))")
    return isempty(parts) ? "selected sites" : join(parts, "; ")
end

function _warn_inconsistent_site_orbitals(group, candidates::Vector{WinOrbital}, requested::Set{String})
    by_site = Dict{String, Set{String}}()
    for orbital in candidates
        push!(get!(by_site, orbital.site_label, Set{String}()), orbital.orbital)
    end
    for site_label in sort!(collect(keys(by_site)))
        missing = setdiff(requested, by_site[site_label])
        isempty(missing) && continue
        @warn(
            "group '$(group.label)': selected site $site_label lacks requested orbitals " *
            "$(sort!(collect(missing))); available = $(sort!(collect(by_site[site_label])))",
        )
    end
    return nothing
end

function build_projection_spec(basis::WinBasis, group_configs)
    groups = ProjectionGroup[]
    for group in group_configs
        (_group_has_site_selector(group) || _group_has_orbital_selector(group)) ||
            error("group '$(group.label)' matched zero orbitals: mode=\"win_groups\" requires at least one selector")
        if !basis.spinors && lowercase(strip(group.spin)) in ("up", "dn", "down")
            error("group '$(group.label)': spin='$(group.spin)' requested but win has spinors=.false.")
        end

        requested = _group_has_orbital_selector(group) ?
            Set(expand_orbitals_and_shells(group.orbitals, group.orbital_shells; label=group.label)) :
            nothing
        candidates = [orb for orb in basis.orbitals if _matches_site(group, orb) && _matches_spin(group, orb)]
        isempty(candidates) && error("group '$(group.label)' matched zero orbitals for $(_selector_description(group))")

        if !isnothing(requested)
            available = Set(orb.orbital for orb in candidates)
            missing = setdiff(requested, available)
            isempty(missing) || error(
                "group '$(group.label)': orbitals $(sort!(collect(missing))) are not defined in win for " *
                "$(_selector_description(group)); available = $(sort!(collect(available)))",
            )
            _warn_inconsistent_site_orbitals(group, candidates, requested)
        end

        selected = WinOrbital[]
        for orbital in candidates
            isnothing(requested) || orbital.orbital in requested || continue
            push!(selected, orbital)
        end
        isempty(selected) && error("group '$(group.label)' matched zero orbitals")
        indices = [orbital.index for orbital in selected]
        atom_count = length(unique(orbital.site_label for orbital in selected))
        push!(groups, ProjectionGroup(group.label, indices, group.color, atom_count))
    end
    return ProjectionSpec(groups, basis.num_wann)
end

end
