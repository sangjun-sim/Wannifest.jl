module Service

using Printf

using ..Bands
using ..Compare
using ..Dos
using ..KPath
using ..LatticeIO
using ..Model: BandResult, BandRunResult, HrData, ReferenceConfig, RunConfig, WsvecTable
using ..ObservablesOam
using ..ObservablesOutput
using ..ObservablesSam
using ..Output
using ..PlotService
using ..Projection
using ..QeXmlIO
using ..ReferenceData
using ..VaspXmlIO
using ..SpinLayout
using ..WannierHrIO
using ..WannierWsvecIO

export run, print_summary

struct Channel
    spin::Symbol
    hr::HrData
    wsvec::Union{Nothing, WsvecTable}
end

function _check_hr_hermiticity(hr; atol::Float64=1e-10)
    maxdiff = WannierHrIO.pair_hermiticity_error(hr)
    if maxdiff > atol
        @warn "hr.dat deviates from Hermitian-pair symmetry by $maxdiff (threshold: $atol)"
        return false
    end
    return true
end

function _read_wsvec(
    path::Union{Nothing, String},
    hr::HrData,
    spin_layout::Symbol,
)::Union{Nothing, WsvecTable}
    isnothing(path) && return nothing
    return WannierWsvecIO.read_wsvec(path; num_wann=hr.num_wann, spin_layout=spin_layout)
end

function _channels(config::RunConfig)
    up_hr = WannierHrIO.read_hr(config.files.hr_path; spin_layout=config.spin.layout)
    up = Channel(:none, up_hr, _read_wsvec(config.files.wsvec_path, up_hr, config.spin.layout))
    isnothing(config.files.hr_down_path) && return (up,)

    down_hr = WannierHrIO.read_hr(config.files.hr_down_path)
    down = Channel(:down, down_hr, _read_wsvec(config.files.wsvec_down_path, down_hr, :qe))
    return (Channel(:up, up.hr, up.wsvec), down)
end

function _compute_band_result(channel::Channel, kpath, config::RunConfig)::BandResult
    projection_spec = Projection.build_projection_spec(config, channel.hr.num_wann)
    basis_transform = Projection.build_basis_rotation_transform(config, channel.hr.num_wann)
    oam_context = ObservablesOam.build_oam_context(
        config,
        channel.hr.num_wann;
        projection_basis_transform=basis_transform,
    )
    sam_context = ObservablesSam.build_sam_context(
        config,
        channel.hr.num_wann;
        projection_basis_transform=basis_transform,
    )
    return Bands.compute_bands(
        channel.hr,
        kpath,
        config;
        wsvec=channel.wsvec,
        projection_spec=projection_spec,
        basis_transform=basis_transform,
        oam_context=oam_context,
        sam_context=sam_context,
    )
end

function print_summary(result::BandRunResult; make_plot::Bool=true, io::IO=stdout)
    SpinLayout.print_summary(result; make_plot=make_plot, io=io)
    comparison = result.comparison
    isnothing(comparison) && return nothing

    println(io, "$(_reference_source_name(result.config.reference.source)) reference-band comparison:")
    comparison_down = result.comparison_down
    if !isnothing(comparison_down)
        @printf(io, "  Up-channel RMS deviation: %.6f eV\n", comparison.overall_rms_eV)
        @printf(io, "  Down-channel RMS deviation: %.6f eV\n", comparison_down.overall_rms_eV)
        @printf(io, "  Overall RMS deviation: %.6f eV\n", Compare.combined_rms(comparison, comparison_down))
        return nothing
    end

    @printf(io, "  Overall RMS deviation: %.6f eV\n", comparison.overall_rms_eV)
    for iband in eachindex(comparison.n_in_window)
        @printf(
            io,
            "  Wannier band %d: n=%d, RMS=%.6f eV, max=%.6f eV\n",
            iband,
            comparison.n_in_window[iband],
            comparison.rms_eV[iband],
            comparison.max_eV[iband],
        )
    end
    return nothing
end

_reference_source_name(source::Symbol) = source == :vasp ? "VASP" : "QE"

function _read_reference_band_data(reference::ReferenceConfig)::ReferenceBandData
    reference.source == :qe_xml &&
        return QeXmlIO.read_qe_band_data(reference.bands_xml_path, reference.scf_xml_path)
    reference.source == :vasp &&
        return VaspXmlIO.read_vasp_band_data(reference.bands_xml_path, reference.scf_xml_path)
    error("unsupported reference.source: $(reference.source)")
end

function run(config::RunConfig; make_plot::Bool=true)::BandRunResult
    channels = _channels(config)
    channel_hermiticity = [_check_hr_hermiticity(channel.hr; atol=config.hermiticity_tol) for channel in channels]
    hermiticity_ok = all(channel_hermiticity)
    lattice = if config.mode in (:bands, :all) && !isempty(config.files.structure_path)
        LatticeIO.read_lattice(config.files.structure_path)
    else
        nothing
    end

    band_result = nothing
    band_result_down = nothing
    compare_band_result = nothing
    compare_band_result_down = nothing
    comparison = nothing
    comparison_down = nothing
    compare_fermi_level = nothing
    if config.mode in (:bands, :all)
        qe_data = nothing
        kpath = if config.reference.enabled
            isnothing(lattice) && error("reference comparison requires run.structure for physical distances")
            qe_data = _read_reference_band_data(config.reference)
            ReferenceData.assert_lattice_match(qe_data, lattice)
            source_name = _reference_source_name(config.reference.source)
            qe_data.lsda && length(channels) == 1 &&
                error("$source_name 참조가 lsda=true(spin-polarized)입니다. run.hr_down으로 up/dn 두 hr 파일을 주십시오")
            !qe_data.lsda && length(channels) == 2 &&
                error("run.hr_down으로 두 채널을 주었지만 $source_name 참조는 lsda=false입니다")
            Compare.qe_kpath(qe_data, lattice, config.reference.k_labels)
        else
            isempty(config.files.kpoints_path) && error("Band mode requires [run] kpoints in input.toml")
            kpd = KPath.parse_kpoints(config.files.kpoints_path; lattice=lattice)
            KPath.generate_kpath(kpd; lattice=lattice)
        end
        channel_results = [_compute_band_result(channel, kpath, config) for channel in channels]
        band_result = channel_results[1]
        band_result_down = length(channel_results) == 2 ? channel_results[2] : nothing
        for (channel, result) in zip(channels, channel_results)
            Output.write_bands_data(
                Output.channel_path(config.output.bands_data, channel.spin),
                result;
                config=config,
                spin=channel.spin,
            )
        end
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
        if config.reference.enabled
            window = isnothing(config.reference.window) ? config.plot.energy_range : config.reference.window
            compare_fermi_level = Compare.reference_fermi_level(qe_data, config.reference)
            if isnothing(band_result_down)
                aligned_wannier, aligned_qe = Compare.aligned_energies(
                    band_result,
                    qe_data,
                    config.reference,
                )
                compare_band_result = BandResult(
                    band_result.kpoints_frac,
                    band_result.distances,
                    aligned_wannier,
                    band_result.tick_positions,
                    band_result.tick_labels,
                    band_result.segment_ranges,
                    band_result.is_physical_distance,
                )
                comparison = Compare.compare_bands(compare_band_result, aligned_qe, window)
                Compare.write_comparison(
                    config.output.compare_data,
                    comparison,
                    window,
                    length(kpath.kpoints),
                )
            else
                aligned_up, aligned_down, qe_up, qe_down = Compare.aligned_energies(
                    band_result,
                    band_result_down,
                    qe_data,
                    config.reference,
                )
                compare_band_result = BandResult(
                    band_result.kpoints_frac,
                    band_result.distances,
                    aligned_up,
                    band_result.tick_positions,
                    band_result.tick_labels,
                    band_result.segment_ranges,
                    band_result.is_physical_distance,
                )
                compare_band_result_down = BandResult(
                    band_result_down.kpoints_frac,
                    band_result_down.distances,
                    aligned_down,
                    band_result_down.tick_positions,
                    band_result_down.tick_labels,
                    band_result_down.segment_ranges,
                    band_result_down.is_physical_distance,
                )
                comparison = Compare.compare_bands(compare_band_result, qe_up, window)
                comparison_down = Compare.compare_bands(compare_band_result_down, qe_down, window)
                Compare.write_comparison(
                    config.output.compare_data,
                    comparison,
                    comparison_down,
                    window,
                    length(kpath.kpoints),
                )
            end
        end
    end

    dos_result = nothing
    if config.mode in (:dos, :all)
        length(channels) == 1 || error("collinear two-channel mode does not support DOS")
        channel = channels[1]
        projection_spec = Projection.build_projection_spec(config, channel.hr.num_wann)
        basis_transform = Projection.build_basis_rotation_transform(config, channel.hr.num_wann)
        dos_result = Dos.run_dos(
            channel.hr,
            config;
            wsvec=channel.wsvec,
            projection_spec=projection_spec,
            basis_transform=basis_transform,
        )
        Output.write_dos_data(config.output.dos_data, dos_result)
        Output.write_pdos_data(config.projection.pdos_data, dos_result)
    end

    PlotService.maybe_plot(
        config,
        band_result,
        band_result_down,
        dos_result,
        compare_band_result,
        compare_band_result_down,
        comparison,
        comparison_down,
        make_plot;
        compare_fermi_level=compare_fermi_level,
    )
    return BandRunResult(
        config,
        band_result,
        band_result_down,
        dos_result,
        hermiticity_ok,
        comparison,
        comparison_down,
    )
end

end
