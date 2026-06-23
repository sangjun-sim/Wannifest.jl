module Service

using ..CenterIO
using ..HrIO
using ..InputIO
using ..Kspace
using ..Model: RunConfig, SuperhamRunResult
using ..PoscarIO: read_poscar
using ..SupercellHam
using ..SupercellGeometry
using ..Validate
using ..WsvecIO
using ..HrFormat: write_hr_blocks_normalized

export run

function _geometry_source(model)::String
    return !isnothing(model.wsvec) && !isnothing(model.centers) ? "wsvec+centers" :
           !isnothing(model.centers) ? String(model.centers.mode) :
           !isnothing(model.wsvec) ? "wsvec-only" :
           "none"
end

function run(
    config::RunConfig;
    kpoint::AbstractVector{<:Real}=Float64[0.0, 0.0, 0.0],
    output_hr::Union{Nothing, AbstractString}=nothing,
)::SuperhamRunResult
    output_path = isnothing(output_hr) ? config.output_hr : String(output_hr)

    structure_cell = read_poscar(config.structure_path)
    lattice = structure_cell.lattice
    model = HrIO.read_hr(config.hr_path; lattice=lattice, spin_layout=config.spin_layout)
    InputIO.validate_manual_orbitals(config, model.num_wann)
    InputIO.warn_if_manual_centers_are_surrogate(config)

    if !isnothing(config.wsvec_path)
        model = WsvecIO.attach_wsvec(
            model,
            WsvecIO.read_wsvec(
                config.wsvec_path;
                num_wann=model.num_wann,
                spin_layout=config.spin_layout,
            ),
        )
    end

    centers = CenterIO.centers_from_config(
        config,
        lattice;
        num_wann=model.num_wann,
        structure_cell=structure_cell,
        spin_layout=config.spin_layout,
    )
    if !isnothing(centers)
        model = CenterIO.attach_centers(model, centers)
    end

    geometry_source = _geometry_source(model)
    primitive_herm = try
        Validate.validate_hermiticity(model; atol=1e-8)
    catch err
        error("Primitive hr.dat hermiticity validation failed: $(sprint(showerror, err))")
    end
    build_result = SupercellHam.build_supercell(
        model,
        SupercellGeometry.from_user_matrix(config.supercell_matrix);
        strict_geometry=config.strict_geometry,
    )
    super_model = build_result.model
    size_result = Validate.validate_supercell_size(model, config.supercell_matrix, super_model)
    fold_diff = Validate.validate_folded_spectrum(model, super_model, config.supercell_matrix; kpoints=[Float64.(kpoint)], atol=1e-8)

    primitive_evals = Kspace.eigenvalues_k(model, kpoint)
    primitive_wsvec_evals = isnothing(model.wsvec) ? nothing : Kspace.eigenvalues_k_ws(model, kpoint)
    supercell_evals = Kspace.eigenvalues_k(super_model, kpoint)

    if !isnothing(output_path)
        write_hr_blocks_normalized(String(output_path), super_model.header, super_model.num_wann, super_model.hoppings)
    end

    return SuperhamRunResult(
        config,
        model,
        super_model,
        output_path,
        geometry_source,
        build_result.report,
        primitive_herm,
        size_result.multiplicity,
        fold_diff,
        primitive_evals,
        primitive_wsvec_evals,
        supercell_evals,
    )
end

end
