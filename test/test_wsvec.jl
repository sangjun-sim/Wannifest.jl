@testset "band strict eigensystem and wsvec policy" begin
    mktempdir() do dir
        nonherm_path = joinpath(dir, "nonhermitian_real_spectrum_hr.dat")
        write(nonherm_path, """
Non-Hermitian real-spectrum model
  2
  1
  1
  0  0  0  1  1   0.000000   0.000000
  0  0  0  2  1   0.000000   0.000000
  0  0  0  1  2   1.000000   0.000000
  0  0  0  2  2   0.000000   0.000000
""")
        nonherm = WannierHrIO.read_hr(nonherm_path)
        nonherm_error = error_message(() -> WannierEigensystem.solve_kpoint(nonherm, [0.0, 0.0, 0.0]))
        @test occursin("not Hermitian", nonherm_error)
    end

    ws_hr = WannierHrIO.read_hr(joinpath(SAMPLE_SPECTRA_DIR, "wsvec_test_hr.dat"))
    full_ws = WannierWsvecIO.read_wsvec(joinpath(SAMPLE_SPECTRA_DIR, "wsvec_test.dat"))
    first_key = first(keys(full_ws.table))
    partial_table = Dict{Tuple{WannierTypes.RKey, Int, Int}, WannierTypes.WsvecEntry}()
    partial_table[first_key] = full_ws.table[first_key]
    partial_ws = WannierTypes.WsvecTable(partial_table)
    missing_error = error_message(
        () -> WannierKspace.hamiltonian_k_wsvec(ws_hr, partial_ws, [0.0, 0.0, 0.0]),
    )
    @test occursin("Missing wsvec entry", missing_error)

    mktempdir() do dir
        ws_path = joinpath(dir, "vasp_block_wsvec.dat")
        write(ws_path, """
VASP block wsvec layout
  0 0 0 3 4
1
  0 0 0
""")
        ws = WannierWsvecIO.read_wsvec(ws_path; num_wann=4, spin_layout=:vasp544)
        @test haskey(ws.table, ((0, 0, 0), 2, 4))
        @test !haskey(ws.table, ((0, 0, 0), 3, 4))

        ws6_path = joinpath(dir, "vasp_block_wsvec_6.dat")
        write(ws6_path, """
VASP block wsvec layout 6
  0 0 0 2 5
1
  0 0 0
""")
        ws6 = WannierWsvecIO.read_wsvec(ws6_path; num_wann=6, spin_layout=:vasp544)
        @test haskey(ws6.table, ((0, 0, 0), 3, 4))
        @test !haskey(ws6.table, ((0, 0, 0), 5, 2))
    end

    avg_hoppings = Dict{WannierTypes.RKey, Matrix{ComplexF64}}(
        (0, 0, 0) => ComplexF64[1.0 + 0.0im;;],
    )
    avg_ndegen = Dict{WannierTypes.RKey, Int}((0, 0, 0) => 7)
    avg_hr = WannierTypes.HrBlocks("wsvec averaging test", 1, 1, avg_hoppings, avg_ndegen, :raw)
    avg_ws = WannierTypes.WsvecTable(Dict{Tuple{WannierTypes.RKey, Int, Int}, WannierTypes.WsvecEntry}(
        ((0, 0, 0), 1, 1) => WannierTypes.WsvecEntry(2, Int[0 2; 0 0; 0 0]),
    ))
    plain_hk = WannierKspace.hamiltonian_k_plain(avg_hr, [0.0, 0.0, 0.0])
    @test plain_hk[1, 1] ≈ (1.0 / 7.0) + 0.0im
    avg_hk = WannierKspace.hamiltonian_k_wsvec(avg_hr, avg_ws, [0.0, 0.0, 0.0])
    @test avg_hk[1, 1] ≈ 1.0 + 0.0im
end

@testset "wsvec generation writer and validation" begin
    lattice = Matrix{Float64}(I, 3, 3)
    centers_one = reshape(Float64[0.0, 0.0, 0.0], 3, 1)
    onsite = Dict{WannierTypes.RKey, Matrix{ComplexF64}}(
        (0, 0, 0) => ComplexF64[1.0 + 0.0im;;],
    )
    onsite_hr = WannierTypes.HrBlocks(
        "trivial wsvec generation",
        1,
        1,
        onsite,
        Dict{WannierTypes.RKey, Int}((0, 0, 0) => 1),
        :raw,
    )

    trivial_ws = WannierWsvecGenerate.generate_wsvec(onsite_hr, lattice, (1, 1, 1), centers_one)
    trivial_entry = trivial_ws.table[((0, 0, 0), 1, 1)]
    @test trivial_entry.n_shift == 1
    @test trivial_entry.shifts == reshape(Int[0, 0, 0], 3, 1)

    boundary_hr = WannierTypes.HrBlocks(
        "boundary wsvec generation",
        1,
        1,
        Dict{WannierTypes.RKey, Matrix{ComplexF64}}((1, 0, 0) => ComplexF64[1.0 + 0.0im;;]),
        Dict{WannierTypes.RKey, Int}((1, 0, 0) => 1),
        :raw,
    )
    boundary_ws = WannierWsvecGenerate.generate_wsvec(boundary_hr, lattice, (2, 1, 1), centers_one; search_size=2)
    boundary_entry = boundary_ws.table[((1, 0, 0), 1, 1)]
    boundary_shifts = Set(
        (boundary_entry.shifts[1, idx], boundary_entry.shifts[2, idx], boundary_entry.shifts[3, idx])
        for idx in 1:boundary_entry.n_shift
    )
    @test boundary_shifts == Set([(0, 0, 0), (-2, 0, 0)])

    boundary_ws_s1 = WannierWsvecGenerate.generate_wsvec(boundary_hr, lattice, (2, 1, 1), centers_one; search_size=1)
    @test boundary_ws_s1.table[((1, 0, 0), 1, 1)].shifts == boundary_entry.shifts

    suspicious_hr = WannierTypes.HrBlocks(
        "suspicious R extent",
        1,
        1,
        Dict{WannierTypes.RKey, Matrix{ComplexF64}}((4, 0, 0) => ComplexF64[1.0 + 0.0im;;]),
        Dict{WannierTypes.RKey, Int}((4, 0, 0) => 1),
        :raw,
    )
    @test_logs (:warn, r"look large") WannierWsvecGenerate.generate_wsvec(
        suspicious_hr,
        lattice,
        (2, 2, 2),
        centers_one,
    )

    y_hr = WannierTypes.HrBlocks(
        "frozen-axis compatibility",
        1,
        1,
        Dict{WannierTypes.RKey, Matrix{ComplexF64}}((0, 1, 0) => ComplexF64[1.0 + 0.0im;;]),
        Dict{WannierTypes.RKey, Int}((0, 1, 0) => 1),
        :raw,
    )
    y_ws = WannierWsvecGenerate.generate_wsvec(y_hr, lattice, (1, 1, 1), centers_one; search_size=1)
    @test y_ws.table[((0, 1, 0), 1, 1)].shifts == reshape(Int[0, -1, 0], 3, 1)
    @test_throws ArgumentError WannierWsvecGenerate.generate_wsvec(
        y_hr,
        lattice,
        (1, 1, 1),
        centers_one;
        frozen_axes=(false, true, false),
    )

    two_center_hr = WannierTypes.HrBlocks(
        "center-sensitive wsvec generation",
        2,
        1,
        Dict{WannierTypes.RKey, Matrix{ComplexF64}}((1, 0, 0) => zeros(ComplexF64, 2, 2)),
        Dict{WannierTypes.RKey, Int}((1, 0, 0) => 1),
        :raw,
    )
    centers_two = Float64[
        0.00 0.75
        0.00 0.00
        0.00 0.00
    ]
    two_center_ws = WannierWsvecGenerate.generate_wsvec(two_center_hr, lattice, (2, 2, 2), centers_two)
    @test two_center_ws.table[((1, 0, 0), 1, 2)].shifts == reshape(Int[-2, 0, 0], 3, 1)
    @test two_center_ws.table[((1, 0, 0), 2, 1)].shifts == reshape(Int[0, 0, 0], 3, 1)

    @test_throws ArgumentError WannierWsvecGenerate.generate_wsvec(
        onsite_hr,
        lattice,
        (1, 1, 1),
        centers_one;
        centers_cart=10.0 .* centers_one .+ 1.0,
    )

    partial_ws = WannierTypes.WsvecTable(Dict{Tuple{WannierTypes.RKey, Int, Int}, WannierTypes.WsvecEntry}(
        ((0, 0, 0), 1, 1) => WannierTypes.WsvecEntry(1, reshape(Int[0, 0, 0], 3, 1)),
    ))
    two_orbital_hr = WannierTypes.HrBlocks(
        "coverage validation",
        2,
        1,
        Dict{WannierTypes.RKey, Matrix{ComplexF64}}((0, 0, 0) => zeros(ComplexF64, 2, 2)),
        Dict{WannierTypes.RKey, Int}((0, 0, 0) => 1),
        :raw,
    )
    coverage = WannierWsvecGenerate.validate_wsvec_coverage(two_orbital_hr, partial_ws)
    @test !coverage.ok
    @test ((0, 0, 0), 1, 2) in coverage.missing_pairs
    @test occursin("wsvec is missing", error_message(() -> WannierWsvecGenerate.assert_wsvec_usable(two_orbital_hr, partial_ws)))

    bad_ws = WannierTypes.WsvecTable(Dict{Tuple{WannierTypes.RKey, Int, Int}, WannierTypes.WsvecEntry}(
        ((0, 0, 0), 1, 1) => WannierTypes.WsvecEntry(2, reshape(Int[0, 0, 0], 3, 1)),
    ))
    @test !WannierWsvecGenerate.validate_wsvec_coverage(onsite_hr, bad_ws).ok
    @test occursin("invalid", error_message(() -> WannierWsvecGenerate.assert_wsvec_usable(onsite_hr, bad_ws)))

    duplicate_shift_ws = WannierTypes.WsvecTable(Dict{Tuple{WannierTypes.RKey, Int, Int}, WannierTypes.WsvecEntry}(
        ((0, 0, 0), 1, 1) => WannierTypes.WsvecEntry(2, Int[0 0; 0 0; 0 0]),
    ))
    duplicate_report = WannierWsvecGenerate.validate_wsvec_coverage(onsite_hr, duplicate_shift_ws)
    @test ((0, 0, 0), 1, 1) in duplicate_report.duplicate_shifts
    @test occursin("duplicate shifts", error_message(() -> WannierWsvecGenerate.assert_wsvec_usable(onsite_hr, duplicate_shift_ws)))

    extra_ws = WannierTypes.WsvecTable(Dict{Tuple{WannierTypes.RKey, Int, Int}, WannierTypes.WsvecEntry}(
        ((0, 0, 0), 1, 1) => WannierTypes.WsvecEntry(1, reshape(Int[0, 0, 0], 3, 1)),
        ((9, 0, 0), 1, 1) => WannierTypes.WsvecEntry(1, reshape(Int[0, 0, 0], 3, 1)),
        ((0, 0, 0), 2, 1) => WannierTypes.WsvecEntry(1, reshape(Int[0, 0, 0], 3, 1)),
    ))
    extra_report = WannierWsvecGenerate.validate_wsvec_coverage(onsite_hr, extra_ws)
    @test ((9, 0, 0), 1, 1) in extra_report.extra_pairs
    @test ((0, 0, 0), 2, 1) in extra_report.out_of_bounds_pairs
    @test occursin("out-of-bounds", error_message(() -> WannierWsvecGenerate.assert_wsvec_usable(onsite_hr, extra_ws)))

    pair_hr = WannierTypes.HrBlocks(
        "hermiticity pair diagnostic",
        1,
        2,
        Dict{WannierTypes.RKey, Matrix{ComplexF64}}(
            (1, 0, 0) => ComplexF64[1.0 + 0.0im;;],
            (-1, 0, 0) => ComplexF64[1.0 + 0.0im;;],
        ),
        Dict{WannierTypes.RKey, Int}((1, 0, 0) => 1, (-1, 0, 0) => 1),
        :raw,
    )
    mismatched_pair_ws = WannierTypes.WsvecTable(Dict{Tuple{WannierTypes.RKey, Int, Int}, WannierTypes.WsvecEntry}(
        ((1, 0, 0), 1, 1) => WannierTypes.WsvecEntry(1, reshape(Int[-2, 0, 0], 3, 1)),
        ((-1, 0, 0), 1, 1) => WannierTypes.WsvecEntry(1, reshape(Int[0, 0, 0], 3, 1)),
    ))
    pair_report = WannierWsvecGenerate.validate_wsvec_coverage(pair_hr, mismatched_pair_ws)
    @test !isempty(pair_report.hermiticity_mismatches)

    one_sided_ws = WannierTypes.WsvecTable(Dict{Tuple{WannierTypes.RKey, Int, Int}, WannierTypes.WsvecEntry}(
        ((1, 0, 0), 1, 1) => WannierTypes.WsvecEntry(1, reshape(Int[0, 0, 0], 3, 1)),
    ))
    one_sided_report = WannierWsvecGenerate.validate_wsvec_coverage(pair_hr, one_sided_ws)
    @test !isempty(one_sided_report.one_sided_pairs)

    ndegen_hr = WannierTypes.HrBlocks(
        "strict ndegen validation",
        1,
        1,
        onsite,
        Dict{WannierTypes.RKey, Int}((0, 0, 0) => 2),
        :raw,
    )
    ndegen_report = WannierWsvecGenerate.assert_wsvec_usable(ndegen_hr, trivial_ws)
    @test (0, 0, 0) in ndegen_report.bad_ndegen
    @test occursin(
        "hr.ndegen == 1",
        error_message(() -> WannierWsvecGenerate.assert_wsvec_usable(ndegen_hr, trivial_ws; require_unit_ndegen=true)),
    )

    sc_hoppings = Dict{SuperhamCore.Model.RKey, Matrix{ComplexF64}}(
        (0, 0, 0) => ComplexF64[1.0 + 0.0im;;],
    )
    sc_ndegen = Dict{SuperhamCore.Model.RKey, Int}((0, 0, 0) => 7)
    sc_ws = SuperhamCore.Model.WsvecTable(Dict{Tuple{SuperhamCore.Model.RKey, Int, Int}, SuperhamCore.Model.WsvecEntry}(
        ((0, 0, 0), 1, 1) => SuperhamCore.Model.WsvecEntry(1, reshape(Int[0, 0, 0], 3, 1)),
    ))
    sc_model = SuperhamCore.Model.HrModel(
        "superham wsvec averaging",
        lattice,
        2π .* inv(lattice)',
        1,
        sc_hoppings,
        sc_ndegen,
        nothing,
        nothing,
    )
    sc_ws_model = SuperhamCore.WsvecIO.attach_wsvec(sc_model, sc_ws)
    @test SuperhamCore.Kspace.hamiltonian_k(sc_ws_model, [0.0, 0.0, 0.0])[1, 1] ≈ (1.0 / 7.0) + 0.0im
    @test SuperhamCore.Kspace.hamiltonian_k_ws(sc_ws_model, [0.0, 0.0, 0.0])[1, 1] ≈ 1.0 + 0.0im

    mktempdir() do dir
        generated_path = joinpath(dir, "nested", "generated_wsvec.dat")
        WannierWsvecGenerate.write_wsvec(
            generated_path,
            boundary_ws;
            mp_grid=(2, 1, 1),
            center_policy=:atomic_assumption,
        )
        roundtrip = WannierWsvecIO.read_wsvec(generated_path)
        @test sort!(collect(keys(roundtrip.table))) == sort!(collect(keys(boundary_ws.table)))
        @test roundtrip.table[((1, 0, 0), 1, 1)].shifts == boundary_entry.shifts
        @test_throws ArgumentError WannierWsvecGenerate.write_wsvec(
            joinpath(dir, "bad_header.dat"),
            boundary_ws;
            header="line one\nline two",
        )

        duplicate_path = joinpath(dir, "duplicate_wsvec.dat")
        write(duplicate_path, """
duplicate wsvec
 0 0 0 1 1
 1
 0 0 0
# comment between entries
 0 0 0 1 1
 1
 0 0 0
""")
        @test occursin("Duplicate wsvec entry", error_message(() -> WannierWsvecIO.read_wsvec(duplicate_path)))
    end
end
