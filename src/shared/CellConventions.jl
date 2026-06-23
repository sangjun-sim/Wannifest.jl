module CellConventions

using LinearAlgebra
using Spglib

using ..CrystalCells: CrystalCell, natoms, wrap_fractional, reciprocal_lattice

export CellConvention, CellTransform, CellFrameBundle
export cell_convention_name, parse_cell_convention
export build_spglib_cell, build_dataset, summarize_dataset
export make_cell_transform, map_kfrac, map_frac_positions
export build_spglib_bundle, spglib_standardized_primitive, spglib_standardized_conventional
export get_cell, get_transform, same_lattice

@enum CellConvention begin
    input_model
    standardized_primitive
    standardized_conventional
    basis_sidecar
end

struct CellTransform
    from_convention::CellConvention
    to_convention::CellConvention
    source_lattice::Matrix{Float64}
    target_lattice::Matrix{Float64}
    k_from_to::Matrix{Float64}
    k_to_from::Matrix{Float64}
    r_from_to::Matrix{Float64}
    r_to_from::Matrix{Float64}
end

struct CellFrameBundle
    input_convention::CellConvention
    input_cell::CrystalCell
    spglib_primitive::Union{Nothing, CrystalCell}
    spglib_conventional::Union{Nothing, CrystalCell}
    transforms::Dict{Tuple{CellConvention, CellConvention}, CellTransform}
    symprec::Float64
end

struct SymmetrySummary
    spacegroup_number::Int
    international_symbol::String
    hall_number::Int
    hall_symbol::String
    choice::String
    pointgroup_symbol::String
end

function cell_convention_name(convention::CellConvention)::String
    return if convention == input_model
        "input_model"
    elseif convention == standardized_primitive
        "standardized_primitive"
    elseif convention == standardized_conventional
        "standardized_conventional"
    elseif convention == basis_sidecar
        "basis_sidecar"
    else
        error("Unsupported CellConvention: $convention")
    end
end

function parse_cell_convention(raw)::CellConvention
    raw isa CellConvention && return raw
    raw isa Symbol && return parse_cell_convention(String(raw))
    raw isa AbstractString || error("cell convention must be a string, symbol, or CellConvention")

    normalized = lowercase(strip(String(raw)))
    normalized = replace(normalized, '-' => '_', ' ' => '_')

    return if normalized == "input_model"
        input_model
    elseif normalized == "standardized_primitive"
        standardized_primitive
    elseif normalized == "standardized_conventional"
        standardized_conventional
    elseif normalized == "basis_sidecar"
        basis_sidecar
    else
        error("Unsupported cell convention: $raw")
    end
end

function build_spglib_cell(cell)
    positions = [Float64.(cell.frac_positions[:, i]) for i in 1:natoms(cell)]
    return Spglib.SpglibCell(Matrix{Float64}(cell.lattice), positions, Int.(cell.species_ids))
end

function build_dataset(
    cell;
    symprec::Float64=1e-5,
    angle_tolerance::Float64=-1.0,
    hall_number::Union{Nothing, Int}=nothing,
)
    angle_tolerance == -1.0 || error("Spglib.jl wrapper currently supports only angle_tolerance = -1.0")
    raw_cell = build_spglib_cell(cell)
    if isnothing(hall_number)
        return Spglib.get_dataset(raw_cell, symprec)
    end
    return Spglib.get_dataset_with_hall_number(raw_cell, hall_number, symprec)
end

function summarize_dataset(dataset)::SymmetrySummary
    return SymmetrySummary(
        Int(dataset.spacegroup_number),
        String(dataset.international_symbol),
        Int(dataset.hall_number),
        String(dataset.hall_symbol),
        String(dataset.choice),
        String(dataset.pointgroup_symbol),
    )
end

function _crystal_cell_from_spglib(raw_cell, template; comment_suffix::String="")
    positions = hcat([Float64.(position) for position in raw_cell.positions]...)
    comment = isempty(comment_suffix) ? String(template.comment) : "$(template.comment) [$comment_suffix]"
    return CrystalCell(
        comment,
        Matrix{Float64}(raw_cell.lattice),
        wrap_fractional(positions),
        copy(template.species_names),
        Int.(raw_cell.atoms),
        String(template.source),
    )
end

function spglib_standardized_conventional(cell; symprec::Float64=1e-5)
    raw = Spglib.standardize_cell(build_spglib_cell(cell), symprec; to_primitive=false, no_idealize=false)
    return _crystal_cell_from_spglib(raw, cell; comment_suffix="standardized conventional")
end

function spglib_standardized_primitive(cell; symprec::Float64=1e-5)
    raw = Spglib.find_primitive(build_spglib_cell(cell), symprec)
    return _crystal_cell_from_spglib(raw, cell; comment_suffix="standardized primitive")
end

function make_cell_transform(
    source_lattice::AbstractMatrix{<:Real},
    target_lattice::AbstractMatrix{<:Real};
    from::CellConvention=input_model,
    to::CellConvention=input_model,
)::CellTransform
    source = Matrix{Float64}(source_lattice)
    target = Matrix{Float64}(target_lattice)
    size(source) == (3, 3) || error("source_lattice must be 3x3")
    size(target) == (3, 3) || error("target_lattice must be 3x3")

    b_source = reciprocal_lattice(source)
    b_target = reciprocal_lattice(target)
    k_from_to = Matrix{Float64}(b_target \ b_source)
    r_from_to = Matrix{Float64}(target \ source)

    return CellTransform(
        from,
        to,
        source,
        target,
        k_from_to,
        inv(k_from_to),
        r_from_to,
        inv(r_from_to),
    )
end

function map_kfrac(transform::CellTransform, k_frac::AbstractVector{<:Real})::Vector{Float64}
    length(k_frac) == 3 || error("k_frac must have length 3")
    return Vector{Float64}(transform.k_from_to * Float64.(k_frac))
end

function map_kfrac(transform::CellTransform, k_frac::NTuple{3, <:Real})::Vector{Float64}
    return map_kfrac(transform, collect(k_frac))
end

function map_kfrac(transform::CellTransform, k_frac::AbstractMatrix{<:Real})::Matrix{Float64}
    size(k_frac, 1) == 3 || error("k_frac matrix must be 3xN")
    return Matrix{Float64}(transform.k_from_to * Float64.(k_frac))
end

function map_frac_positions(transform::CellTransform, xs::AbstractVector{<:Real})::Vector{Float64}
    length(xs) == 3 || error("fractional position must have length 3")
    return Vector{Float64}(transform.r_from_to * Float64.(xs))
end

function map_frac_positions(transform::CellTransform, xs::NTuple{3, <:Real})::Vector{Float64}
    return map_frac_positions(transform, collect(xs))
end

function map_frac_positions(transform::CellTransform, xs::AbstractMatrix{<:Real})::Matrix{Float64}
    size(xs, 1) == 3 || error("fractional position matrix must be 3xN")
    return Matrix{Float64}(transform.r_from_to * Float64.(xs))
end

function same_lattice(a::AbstractMatrix{<:Real}, b::AbstractMatrix{<:Real}; atol::Float64=1e-8, rtol::Float64=1e-8)::Bool
    size(a) == (3, 3) || error("lattice a must be 3x3")
    size(b) == (3, 3) || error("lattice b must be 3x3")
    return isapprox(Matrix{Float64}(a), Matrix{Float64}(b); atol=atol, rtol=rtol)
end

function build_spglib_bundle(
    cell;
    symprec::Float64=1e-5,
    input_convention::CellConvention=input_model,
)::CellFrameBundle
    input_cell = CrystalCell(
        String(cell.comment),
        Matrix{Float64}(cell.lattice),
        Matrix{Float64}(cell.frac_positions),
        copy(cell.species_names),
        Int.(cell.species_ids),
        String(cell.source),
    )

    primitive = spglib_standardized_primitive(input_cell; symprec=symprec)
    conventional = spglib_standardized_conventional(input_cell; symprec=symprec)

    cells = Dict{CellConvention, CrystalCell}(input_convention => input_cell)
    cells[standardized_primitive] = primitive
    cells[standardized_conventional] = conventional

    transforms = Dict{Tuple{CellConvention, CellConvention}, CellTransform}()
    for (from_convention, source_cell) in cells
        for (to_convention, target_cell) in cells
            transforms[(from_convention, to_convention)] = make_cell_transform(
                source_cell.lattice,
                target_cell.lattice;
                from=from_convention,
                to=to_convention,
            )
        end
    end

    return CellFrameBundle(
        input_convention,
        input_cell,
        primitive,
        conventional,
        transforms,
        symprec,
    )
end

function get_cell(bundle::CellFrameBundle, convention::CellConvention)
    if convention == bundle.input_convention
        return bundle.input_cell
    elseif convention == standardized_primitive
        isnothing(bundle.spglib_primitive) && error("CellFrameBundle has no standardized_primitive cell")
        return bundle.spglib_primitive
    elseif convention == standardized_conventional
        isnothing(bundle.spglib_conventional) && error("CellFrameBundle has no standardized_conventional cell")
        return bundle.spglib_conventional
    end
    error("CellFrameBundle does not store convention $(cell_convention_name(convention))")
end

get_cell(bundle::CellFrameBundle, raw) = get_cell(bundle, parse_cell_convention(raw))

function get_transform(bundle::CellFrameBundle, from, to)::CellTransform
    from_conv = parse_cell_convention(from)
    to_conv = parse_cell_convention(to)
    key = (from_conv, to_conv)
    haskey(bundle.transforms, key) || error(
        "CellFrameBundle has no transform $(cell_convention_name(from_conv)) -> $(cell_convention_name(to_conv))",
    )
    return bundle.transforms[key]
end

end
