module Service

using ..BasisSource
using ..BondGeometry
using ..Diagnostics
using ..FluxCompiler
using ..HrFormat
using ..LatticeIO
using ..Model: FluxConfig, FluxRunResult
using ..PairChecks
using ..PlotFlux3D
using ..WannierHrIO

export default_output_hr, default_html_path, default_diagnostic_path, run

function _stem(path::AbstractString)
    name = basename(path)
    for suffix in ("_hr.dat", ".dat")
        endswith(name, suffix) && return name[1:end - length(suffix)]
    end
    return splitext(name)[1]
end

function default_output_hr(config::FluxConfig)
    dir = dirname(config.files.hr_path)
    return joinpath(dir, string(_stem(config.files.hr_path), "_flux_hr.dat"))
end

function default_html_path(config::FluxConfig)
    return joinpath(dirname(config.files.hr_path), "outputs", "plots", string(_stem(config.files.hr_path), "_flux.html"))
end

function default_diagnostic_path(config::FluxConfig)
    return joinpath(dirname(config.files.hr_path), "outputs", "diagnostic", string(_stem(config.files.hr_path), "_flux_diagnostic.tsv"))
end

function run(
    config::FluxConfig;
    output_hr::Union{Nothing, AbstractString}=nothing,
    html_path::Union{Nothing, AbstractString}=nothing,
    diagnostic_path::Union{Nothing, AbstractString}=nothing,
    make_html::Bool=config.plot.interactive,
    make_diagnostic::Bool=config.diagnostic.enabled,
    validate_roundtrip::Bool=false,
)::FluxRunResult
    hr = WannierHrIO.read_hr(config.files.hr_path; spin_layout=config.spin.layout)
    basis, lattice = if !isnothing(config.files.win_path)
        entries, _ = BasisSource.read_win_flux_basis(config.files.win_path; spin_layout=config.spin.layout)
        entries, LatticeIO.read_wannier_win(config.files.win_path).real_lattice
    else
        poscar_path = config.files.poscar_path
        isnothing(poscar_path) && error("flux.run requires either win or poscar")
        entries, cell = BasisSource.read_poscar_flux_basis(
            poscar_path,
            hr.num_wann;
            orbitals_per_atom=config.basis.orbitals_per_atom,
            orbitals_per_species_group=config.basis.orbitals_per_species_group,
        )
        entries, cell.lattice
    end
    length(basis) == hr.num_wann ||
        error("basis source has $(length(basis)) orbitals but hr has num_wann=$(hr.num_wann)")
    pairs = BondGeometry.enumerate_pairs(
        basis,
        lattice;
        search_bounds=config.geometry.search_bounds,
        distance_tol=config.geometry.distance_tol,
    )
    hops = WannierHrIO.normalized_hoppings(hr)
    edges = FluxCompiler.apply_flux_terms!(hops, basis, pairs, lattice, config, hr.num_wann)
    FluxCompiler.complete_hermiticity!(hops)

    diagnostic = if make_diagnostic
        Diagnostics.run_diagnostics(edges, basis, config.diagnostic)
    else
        nothing
    end

    out_hr = String(isnothing(output_hr) ? default_output_hr(config) : output_hr)
    HrFormat.write_hr_blocks_normalized(out_hr, string(hr.header, " + flux terms"), hr.num_wann, hops)

    html_out = if make_html
        out_html = String(isnothing(html_path) ? default_html_path(config) : html_path)
        plot_bounds = isnothing(config.plot.cell_bounds) ? config.geometry.search_bounds : config.plot.cell_bounds
        PlotFlux3D.write_flux_html(out_html, lattice, basis, edges, plot_bounds, config.plot)
    else
        nothing
    end

    diagnostic_out = if !isnothing(diagnostic)
        out_diagnostic = String(isnothing(diagnostic_path) ? default_diagnostic_path(config) : diagnostic_path)
        Diagnostics.write_diagnostic_tsv(out_diagnostic, diagnostic)
    else
        nothing
    end

    if validate_roundtrip
        parsed = WannierHrIO.read_hr(out_hr; spin_layout=config.spin.layout)
        PairChecks.check_hr_pair_symmetry(WannierHrIO.normalized_hoppings(parsed))
    end

    return FluxRunResult(config, out_hr, html_out, diagnostic_out, edges, diagnostic, validate_roundtrip)
end

end
