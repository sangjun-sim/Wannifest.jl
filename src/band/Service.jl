module Service

using ..Bands
using ..Dos
using ..KPath
using ..LatticeIO
using ..Model: BandRunResult, RunConfig
using ..ObservablesOam
using ..ObservablesOutput
using ..ObservablesSam
using ..Output
using ..PlotService
using ..Projection
using ..SpinLayout
using ..WannierHrIO
using ..WannierWsvecGenerate
using ..WannierWsvecIO

export run, print_summary

_validate_wsvec(hr, wsvec) = WannierWsvecGenerate.assert_wsvec_usable(hr, wsvec)

function _check_hr_hermiticity(hr; atol::Float64=1e-10)
    maxdiff = WannierHrIO.pair_hermiticity_error(hr)
    if maxdiff > atol
        @warn "hr.dat deviates from Hermitian-pair symmetry by $maxdiff (threshold: $atol)"
        return false
    end
    return true
end

function print_summary(result::BandRunResult; make_plot::Bool=true, io::IO=stdout)
    return SpinLayout.print_summary(result; make_plot=make_plot, io=io)
end

function run(config::RunConfig; make_plot::Bool=true)::BandRunResult
    hr = WannierHrIO.read_hr(config.files.hr_path; spin_layout=config.spin.layout)
    hermiticity_ok = _check_hr_hermiticity(hr; atol=config.hermiticity_tol)
    projection_spec = Projection.build_projection_spec(config, hr.num_wann)
    basis_transform = Projection.build_basis_rotation_transform(config, hr.num_wann)
    oam_context = ObservablesOam.build_oam_context(
        config,
        hr.num_wann;
        projection_basis_transform=basis_transform,
    )
    sam_context = ObservablesSam.build_sam_context(
        config,
        hr.num_wann;
        projection_basis_transform=basis_transform,
    )

    wsvec_path = config.files.wsvec_path
    wsvec = isnothing(wsvec_path) ? nothing : WannierWsvecIO.read_wsvec(
        wsvec_path;
        num_wann=hr.num_wann,
        spin_layout=config.spin.layout,
    )
    isnothing(wsvec) || _validate_wsvec(hr, wsvec)
    lattice = if config.mode in (:bands, :all) && !isempty(config.files.structure_path)
        LatticeIO.read_lattice(config.files.structure_path)
    else
        nothing
    end

    band_result = nothing
    if config.mode in (:bands, :all)
        isempty(config.files.kpoints_path) && error("Band mode requires [run] kpoints in input.toml")
        kpd = KPath.parse_kpoints(config.files.kpoints_path; lattice=lattice)
        kpath = KPath.generate_kpath(kpd; lattice=lattice)
        band_result = Bands.compute_bands(
            hr,
            kpath,
            config;
            wsvec=wsvec,
            projection_spec=projection_spec,
            basis_transform=basis_transform,
            oam_context=oam_context,
            sam_context=sam_context,
        )
        Output.write_bands_data(config.output.bands_data, band_result; config=config)
        Output.write_projection_weights_data(config.projection.weights_data, band_result)
        ObservablesOutput.write_oam_data(
            ObservablesOutput.oam_data_path(config.output.bands_data),
            band_result,
            config,
        )
        ObservablesOutput.write_sam_data(
            ObservablesOutput.sam_data_path(config.output.bands_data),
            band_result,
            config,
        )
    end

    dos_result = nothing
    if config.mode in (:dos, :all)
        dos_result = Dos.run_dos(
            hr,
            config;
            wsvec=wsvec,
            projection_spec=projection_spec,
            basis_transform=basis_transform,
        )
        Output.write_dos_data(config.output.dos_data, dos_result)
        Output.write_pdos_data(config.projection.pdos_data, dos_result)
    end

    PlotService.maybe_plot(config, band_result, dos_result, make_plot)
    return BandRunResult(config, band_result, dos_result, hermiticity_ok)
end

end
