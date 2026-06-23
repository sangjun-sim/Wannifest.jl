module CenterIO

using LinearAlgebra

using ..CrystalCells: wrap_fractional
using ..Model: CenterLoadSpec, CenterTable, HrModel, OrbitalSpec, RunConfig
using ..PoscarIO: read_poscar
using ..SpinLayout

export read_centres_xyz, build_manual_centers, build_atomic_centers, canonicalize_centers
export center_load_spec, load_centers, centers_from_config, attach_centers

function read_centres_xyz(path::AbstractString; num_wann::Int, lattice::Matrix{Float64})
    open(path, "r") do io
        total = parse(Int, strip(readline(io)))
        total >= num_wann || error("centres.xyz contains only $total entries, expected at least $num_wann")
        readline(io)

        centers_cart = Matrix{Float64}(undef, 3, num_wann)
        centers_frac = Matrix{Float64}(undef, 3, num_wann)
        labels = Vector{String}(undef, num_wann)

        for j in 1:num_wann
            fields = split(strip(readline(io)))
            length(fields) >= 4 || error("Malformed centres.xyz line for center $j")
            labels[j] = String(fields[1])
            cart = parse.(Float64, fields[2:4])
            frac = wrap_fractional(lattice \ cart)
            centers_frac[:, j] = frac
            centers_cart[:, j] = lattice * frac
        end

        return CenterTable(centers_frac, centers_cart, labels, abspath(path), :wannier)
    end
end

function _center_table_from_specs(
    specs::Vector{OrbitalSpec},
    lattice::Matrix{Float64},
    source::AbstractString,
    mode::Symbol,
)
    n = length(specs)
    centers_frac = Matrix{Float64}(undef, 3, n)
    centers_cart = Matrix{Float64}(undef, 3, n)
    labels = Vector{String}(undef, n)

    for (j, spec) in enumerate(specs)
        frac = wrap_fractional(Float64[spec.center_frac[1], spec.center_frac[2], spec.center_frac[3]])
        centers_frac[:, j] = frac
        centers_cart[:, j] = lattice * frac
        labels[j] = spec.label
    end

    return CenterTable(centers_frac, centers_cart, labels, String(source), mode)
end

function build_manual_centers(
    lattice::Matrix{Float64},
    specs::Vector{OrbitalSpec};
    source::AbstractString="input.toml",
)
    return _center_table_from_specs(specs, lattice, source, :manual_centers)
end

function build_atomic_centers(
    poscar_path::AbstractString,
    lattice::Matrix{Float64},
    specs::Vector{OrbitalSpec};
    source::AbstractString="input.toml",
    cell=nothing,
)
    @warn (
        "geometry mode atomic_assumption uses an atomic-centered approximation: " *
        "orbital centers are taken from atom positions in the structure file " *
        "(POSCAR or wannier90.win-derived atoms) " *
        "rather than from converged Wannier centers. wsvec/MDRS distances and " *
        "geometry-dependent observables are approximate unless the Wannier centers " *
        "are actually atom-centered."
    ) structure=abspath(poscar_path) source=source
    cell = isnothing(cell) ? read_poscar(poscar_path) : cell
    natoms = size(cell.frac_positions, 2)
    rewritten = OrbitalSpec[]
    for spec in specs
        isnothing(spec.atom_index) && error("atomic_assumption requires atom_index for every orbital")
        1 <= spec.atom_index <= natoms || error("atom_index=$(spec.atom_index) out of range 1:$natoms")
        pos = cell.frac_positions[:, spec.atom_index]
        push!(rewritten, OrbitalSpec(spec.label, (pos[1], pos[2], pos[3]), spec.atom_index))
    end
    return _center_table_from_specs(rewritten, lattice, source, :atomic_assumption)
end

function canonicalize_centers(centers::CenterTable, spin_layout)::CenterTable
    num_wann = size(centers.centers_frac, 2)
    index_map = SpinLayout.source_to_canonical_indices(num_wann, spin_layout)
    index_map == collect(1:num_wann) && return centers

    centers_frac = similar(centers.centers_frac)
    centers_cart = similar(centers.centers_cart)
    labels = similar(centers.labels)
    for source in 1:num_wann
        target = index_map[source]
        centers_frac[:, target] = centers.centers_frac[:, source]
        centers_cart[:, target] = centers.centers_cart[:, source]
        labels[target] = centers.labels[source]
    end
    return CenterTable(centers_frac, centers_cart, labels, centers.source, centers.mode)
end

function center_load_spec(cfg::RunConfig)::CenterLoadSpec
    return CenterLoadSpec(
        cfg.geometry_mode,
        cfg.centres_path,
        cfg.manual_num_wann,
        cfg.orbital_specs,
        cfg.structure_path,
        "superham",
    )
end

function load_centers(
    spec::CenterLoadSpec,
    lattice::Matrix{Float64};
    num_wann::Union{Nothing, Int}=nothing,
    structure_cell=nothing,
    spin_layout,
)
    centers = nothing
    if !isnothing(spec.centres_path)
        target_num_wann = isnothing(spec.manual_num_wann) ? num_wann : spec.manual_num_wann
        if isnothing(target_num_wann)
            target_num_wann = length(spec.orbital_specs)
        end
        target_num_wann > 0 || error("centres_path was provided but num_wann could not be inferred")
        centers = read_centres_xyz(spec.centres_path; num_wann=target_num_wann, lattice=lattice)
    elseif spec.mode == :manual_centers
        isempty(spec.orbital_specs) &&
            error("geometry mode :manual_centers requires geometry.orbitals in [$(spec.source_context).geometry]")
        centers = build_manual_centers(lattice, spec.orbital_specs; source=spec.structure_path)
    elseif spec.mode == :atomic_assumption
        isempty(spec.orbital_specs) &&
            error("geometry mode :atomic_assumption requires geometry.orbitals in [$(spec.source_context).geometry]")
        centers = build_atomic_centers(
            spec.structure_path,
            lattice,
            spec.orbital_specs;
            source=spec.structure_path,
            cell=structure_cell,
        )
    end
    return isnothing(centers) ? nothing : canonicalize_centers(centers, spin_layout)
end

function centers_from_config(
    cfg::RunConfig,
    lattice::Matrix{Float64};
    num_wann::Union{Nothing, Int}=nothing,
    structure_cell=nothing,
    spin_layout=cfg.spin_layout,
)
    return load_centers(
        center_load_spec(cfg),
        lattice;
        num_wann=num_wann,
        structure_cell=structure_cell,
        spin_layout=spin_layout,
    )
end

function attach_centers(model::HrModel, centers::CenterTable)
    return HrModel(
        model.header,
        model.lattice,
        model.reciprocal,
        model.num_wann,
        model.hoppings,
        model.ndegen,
        model.wsvec,
        centers,
        normalization=model.blocks.normalization,
    )
end

end
