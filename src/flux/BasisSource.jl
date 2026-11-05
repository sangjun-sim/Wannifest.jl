module BasisSource

using ..Model: FluxBasisEntry
using ..PoscarIO
using ..SpinLayout
using ..Win90Basis

export read_poscar_flux_basis, read_win_flux_basis, entry_label

function read_win_flux_basis(path::AbstractString; spin_layout=SpinLayout.DEFAULT_LAYOUT)
    basis = Win90Basis.read_win_basis(path; spin_layout=spin_layout)
    entries = FluxBasisEntry[
        FluxBasisEntry(
            orbital.index,
            orbital.site_label,
            orbital.species,
            orbital.orbital,
            orbital.spin,
            orbital.center_frac,
        ) for orbital in basis.orbitals
    ]
    return entries, basis
end

function _poscar_atom_labels(cell)
    counters = Dict{String, Int}()
    labels = String[]
    for species in PoscarIO.species_labels(cell)
        counters[species] = get(counters, species, 0) + 1
        sep = !isempty(species) && isdigit(last(species)) ? "_" : ""
        push!(labels, string(species, sep, counters[species]))
    end
    return labels
end

function _count_for_species(counts::Dict{String, Int}, species::String)
    if haskey(counts, species)
        return counts[species]
    elseif haskey(counts, "*")
        return counts["*"]
    end
    error("POSCAR fallback needs flux.basis.orbitals_per_atom entry for species '$species'")
end

function _poscar_orbital_counts(
    cell,
    num_wann::Int,
    counts::Dict{String, Int},
    group_counts::Vector{Int},
)
    species = PoscarIO.species_labels(cell)
    if !isempty(group_counts)
        length(group_counts) == length(cell.species_names) ||
            error("flux.basis.orbitals_per_atom group count length must match POSCAR species groups")
        per_atom = [group_counts[cell.species_ids[atom_index]] for atom_index in 1:PoscarIO.natoms(cell)]
        sum(per_atom) == num_wann ||
            error("flux.basis.orbitals_per_atom expands to $(sum(per_atom)) orbitals, but hr num_wann=$num_wann")
        return per_atom
    end
    if isempty(counts)
        nat = PoscarIO.natoms(cell)
        num_wann % nat == 0 ||
            error("POSCAR fallback cannot infer orbitals per atom: num_wann=$num_wann, natoms=$nat")
        return fill(div(num_wann, nat), nat)
    end
    per_atom = [_count_for_species(counts, item) for item in species]
    sum(per_atom) == num_wann ||
        error("flux.basis.orbitals_per_atom expands to $(sum(per_atom)) orbitals, but hr num_wann=$num_wann")
    return per_atom
end

function read_poscar_flux_basis(
    path::AbstractString,
    num_wann::Int;
    orbitals_per_atom::Dict{String, Int}=Dict{String, Int}(),
    orbitals_per_species_group::Vector{Int}=Int[],
)
    cell = PoscarIO.read_poscar(path)
    species = PoscarIO.species_labels(cell)
    site_labels = _poscar_atom_labels(cell)
    counts = _poscar_orbital_counts(cell, num_wann, orbitals_per_atom, orbitals_per_species_group)
    entries = FluxBasisEntry[]
    for atom_index in 1:PoscarIO.natoms(cell)
        center = (
            cell.frac_positions[1, atom_index],
            cell.frac_positions[2, atom_index],
            cell.frac_positions[3, atom_index],
        )
        for local_index in 1:counts[atom_index]
            push!(entries, FluxBasisEntry(
                length(entries) + 1,
                site_labels[atom_index],
                species[atom_index],
                string("orb", local_index),
                :unpolarized,
                center,
            ))
        end
    end
    return entries, cell
end

function entry_label(entry::FluxBasisEntry)
    spin = entry.spin == :unpolarized ? "" : string(":", entry.spin)
    return string(entry.site_label, ":", entry.orbital, spin)
end

end
