@testset "orbital projection core and selectors" begin
    evecs = Matrix{ComplexF64}(I, 4, 4)
    spec = OrbitalProjection.ProjectionSpec(
        [OrbitalProjection.ProjectionGroup("A", [1, 2], "blue")],
        4,
    )
    weights = OrbitalProjection.weights_for_eigenvectors(spec, evecs)
    @test weights[:, 1] ≈ [1.0, 1.0, 0.0, 0.0]
    @test spec.disjoint
    @test !spec.covers_all

    partition = OrbitalProjection.ProjectionSpec(
        [
            OrbitalProjection.ProjectionGroup("A", [1, 2], "blue"),
            OrbitalProjection.ProjectionGroup("B", [3, 4], "red"),
        ],
        4,
    )
    @test partition.disjoint
    @test partition.covers_all
    @test all(sum(OrbitalProjection.weights_for_eigenvectors(partition, evecs); dims=2) .≈ 1.0)
    unitary = Matrix(qr(randn(ComplexF64, 4, 4)).Q)
    @test all(sum(OrbitalProjection.weights_for_eigenvectors(partition, unitary); dims=2) .≈ 1.0)

    @test_throws ErrorException OrbitalProjection.ProjectionSpec(OrbitalProjection.ProjectionGroup[], 4)
    @test_throws ErrorException OrbitalProjection.ProjectionSpec(
        [OrbitalProjection.ProjectionGroup("A", Int[], "blue")],
        4,
    )
    @test_throws ErrorException OrbitalProjection.ProjectionSpec(
        [OrbitalProjection.ProjectionGroup("A", [0, 1], "blue")],
        4,
    )
    @test_throws ErrorException OrbitalProjection.ProjectionSpec(
        [OrbitalProjection.ProjectionGroup("A", [1, 5], "blue")],
        4,
    )
    @test_throws ErrorException OrbitalProjection.ProjectionSpec(
        [OrbitalProjection.ProjectionGroup("A", [1, 1], "blue")],
        4,
    )
    @test_throws ErrorException OrbitalProjection.ProjectionSpec(
        [
            OrbitalProjection.ProjectionGroup("A", [1], "blue"),
            OrbitalProjection.ProjectionGroup("A", [2], "red"),
        ],
        4,
    )
    @test_throws ErrorException OrbitalProjection.weights_for_eigenvectors(spec, randn(ComplexF64, 3, 3))
    @test OrbitalProjection.ProjectionGroup("A", [2, 1], "blue").atom_count == 1
    @test OrbitalProjection.ProjectionGroup("A", [1, 2], "blue", 2).atom_count == 2
    @test_throws ErrorException OrbitalProjection.ProjectionGroup("A", [1, 2], "blue", 0)

    base_cfg = InputIO.read_input(joinpath(GRAPHENE_EXAMPLE_DIR, "input.band.toml"))
    @test Model.ProjectionConfig === BandCore.ProjectionModel.ProjectionConfig
    @test Model.ProjectionGroupConfig === BandCore.ProjectionModel.ProjectionGroupConfig
    @test Model.ProjectedDosResult === BandCore.ProjectionModel.ProjectedDosResult

    mktempdir() do dir
        hr_path = joinpath(dir, "empty_hr.dat")
        win_path = joinpath(dir, "mo_basis.win")
        input_path = joinpath(dir, "compact_groups.toml")
        write(hr_path, "")
        write(win_path, """
num_wann = 9

begin unit_cell_cart
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
end unit_cell_cart

begin atoms_frac
Mo 0.0 0.0 0.0
end atoms_frac

begin projections
Mo : s; p; d
end projections
""")
        write(input_path, """
[band.run]
mode = "dos"
hr = "empty_hr.dat"
verbose = false

[band.energy]
shift = 0.0

[band.projection]
enabled = true
mode = "win_groups"
win = "mo_basis.win"
groups = [
  ["Mo_sp", ["Mo"], ["s", "p"], "red"],
  ["Mo_t2g", ["Mo"], ["t2g"], "blue"],
]
""")
        compact_cfg = InputIO.read_input(input_path)
        @test compact_cfg.projection.color_group == ["Mo_sp", "Mo_t2g"]
        @test compact_cfg.projection.weights_data ==
            joinpath(dir, "outputs", "data", "projection", "band_projection_weights.dat")
        @test compact_cfg.projection.pdos_data ==
            joinpath(dir, "outputs", "data", "projection", "pdos.dat")
        @test compact_cfg.projection.projected_bands_plot ==
            joinpath(dir, "outputs", "plots", "projection", "bands_projected.png")
        @test compact_cfg.projection.pdos_plot ==
            joinpath(dir, "outputs", "plots", "projection", "pdos.png")
        @test compact_cfg.projection.projected_combined_plot ==
            joinpath(dir, "outputs", "plots", "projection", "band_pdos.png")
        @test [group.label for group in compact_cfg.projection.groups] == ["Mo_sp", "Mo_t2g"]
        @test compact_cfg.projection.groups[1].species == ["Mo"]
        @test compact_cfg.projection.groups[1].orbitals == ["s", "p"]
        @test compact_cfg.projection.groups[1].color == "red"
        @test compact_cfg.projection.groups[2].orbitals == ["t2g"]
        @test compact_cfg.projection.groups[2].color == "blue"
        compact_spec = Win90Basis.build_projection_spec(
            Win90Basis.read_win_basis(win_path),
            compact_cfg.projection.groups,
        )
        @test [group.label for group in compact_spec.groups] == ["Mo_sp", "Mo_t2g"]
        @test length(compact_spec.groups[1].indices) == 4
        @test length(compact_spec.groups[2].indices) == 3
    end

    mktempdir() do dir
        win_path = joinpath(dir, "wrapped.win")
        write(win_path, """
num_wann = 1

begin unit_cell_cart
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
end unit_cell_cart

begin atoms_frac
C -1.0e-15 0.999999999999999 1.000000000000001
end atoms_frac

begin projections
C : s
end projections
""")
        source = BandCore.Win90ProjectionIO.read_win_projection_source(win_path)
        @test source.species_atoms["C"][1] == (0.0, 0.0, 0.0)
        @test source.seeds[1].center_frac == (0.0, 0.0, 0.0)
    end

    index_projection = Model.ProjectionConfig(
        true,
        :index_groups,
        ["A"],
        :colorbar,
        "viridis",
        String[],
        9.0,
        1.0,
        "weights.dat",
        "projected_bands.png",
        "pdos.dat",
        "pdos.png",
        "combined.png",
        nothing,
        [
            Model.ProjectionGroupConfig("A", "blue", [1, 2], String[], Int[], String[], String[], String[], "any"),
        ],
    )
    index_cfg = Model.RunConfig(
        base_cfg.mode,
        base_cfg.files,
        base_cfg.output,
        base_cfg.dos,
        base_cfg.energy,
        base_cfg.plot,
        base_cfg.band_plot,
        base_cfg.dos_plot,
        base_cfg.combined_plot,
        base_cfg.spin,
        index_projection,
        base_cfg.hermiticity_tol,
        base_cfg.verbose,
    )
    projection_spec = Projection.build_projection_spec(index_cfg, 4)
    @test index_cfg.projection.color_group == ["A"]
    @test projection_spec.groups[1].indices == [1, 2]
    @test projection_spec.groups[1].atom_count == 1
    @test Projection.band_projection_result(projection_spec, zeros(Float64, 1, 4, 1)).labels == ["A"]
    rotated_weights = Projection.projection_weights(
        OrbitalProjection.ProjectionSpec(
            [OrbitalProjection.ProjectionGroup("px", [2], "blue")],
            3,
        ),
        Matrix{ComplexF64}(I, 3, 3);
        basis_transform=ComplexF64[
            1.0 0.0 0.0
            0.0 0.0 1.0
            0.0 -1.0 0.0
        ],
    )
    @test rotated_weights[:, 1] ≈ [0.0, 0.0, 1.0]
    per_atom_spec = OrbitalProjection.ProjectionSpec(
        [
            OrbitalProjection.ProjectionGroup("A", [1, 2], "blue", 2),
            OrbitalProjection.ProjectionGroup("B", [3], "red", 1),
        ],
        4,
    )
    per_atom_pdos = Projection.projected_dos_result(per_atom_spec, [4.0 3.0; 8.0 6.0])
    @test per_atom_pdos.atom_counts == [2, 1]
    @test per_atom_pdos.pdos == [2.0 3.0; 4.0 6.0]

    duplicate_projection = Model.ProjectionConfig(
        true,
        :index_groups,
        ["A"],
        :colorbar,
        "viridis",
        String[],
        9.0,
        1.0,
        "weights.dat",
        "projected_bands.png",
        "pdos.dat",
        "pdos.png",
        "combined.png",
        nothing,
        [
            Model.ProjectionGroupConfig("A", "blue", [1], String[], Int[], String[], String[], String[], "any"),
            Model.ProjectionGroupConfig("A", "red", [2], String[], Int[], String[], String[], String[], "any"),
        ],
    )
    duplicate_cfg = Model.RunConfig(
        base_cfg.mode,
        base_cfg.files,
        base_cfg.output,
        base_cfg.dos,
        base_cfg.energy,
        base_cfg.plot,
        base_cfg.band_plot,
        base_cfg.dos_plot,
        base_cfg.combined_plot,
        base_cfg.spin,
        duplicate_projection,
        base_cfg.hermiticity_tol,
        base_cfg.verbose,
    )
    @test occursin("duplicate group labels", error_message(() -> Projection.build_projection_spec(duplicate_cfg, 4)))

    mktempdir() do dir
        write(joinpath(dir, "empty_hr.dat"), "")
        input_path = joinpath(dir, "removed_orbital_groups.toml")
        write(input_path, """
[run]
mode = "dos"
hr = "empty_hr.dat"

[projection]
enabled = true
mode = "orbital_groups"
groups = [["Ru_s", ["Ru"], ["s"], "blue"]]
""")
        @test occursin(
            "Unsupported projection.mode",
            error_message(() -> InputIO.read_input(input_path)),
        )

        basis_key_input = joinpath(dir, "removed_basis_key.toml")
        write(basis_key_input, """
[run]
mode = "dos"
hr = "empty_hr.dat"

[projection]
enabled = true
mode = "index_groups"
basis = "basis.toml"
groups = [["A", [1], "blue"]]
""")
        @test occursin(
            "Unsupported projection option(s): basis",
            error_message(() -> InputIO.read_input(basis_key_input)),
        )
    end

    ru_basis = Win90Basis.read_win_basis(joinpath(RUCL3_EXAMPLE_DIR, "wannier90.win"))
    @test ru_basis.num_wann == 6
    @test [orb.orbital for orb in ru_basis.orbitals[1:3]] == ["dxz", "dyz", "dxy"]
    ru_t2g = Model.ProjectionGroupConfig(
        "Ru1_t2g",
        "#1f77b4",
        Int[],
        ["Ru"],
        [1],
        String[],
        ["dxz", "dxy", "dyz"],
        String[],
        "any",
    )
    ru_spec = Win90Basis.build_projection_spec(ru_basis, [ru_t2g])
    @test ru_spec.groups[1].indices == [1, 2, 3]
    @test ru_spec.groups[1].atom_count == 1

    ru_p = Model.ProjectionGroupConfig(
        "Ru_p",
        "#1f77b4",
        Int[],
        ["Ru"],
        Int[],
        String[],
        String[],
        ["p"],
        "any",
    )
    @test occursin("not defined in win", error_message(() -> Win90Basis.build_projection_spec(ru_basis, [ru_p])))

    la_f = Model.ProjectionGroupConfig(
        "La_f",
        "#1f77b4",
        Int[],
        ["La"],
        Int[],
        String[],
        String[],
        ["f"],
        "any",
    )
    @test occursin("matched zero orbitals", error_message(() -> Win90Basis.build_projection_spec(ru_basis, [la_f])))

    mos2_basis = Win90Basis.read_win_basis(joinpath(MOS2_EXAMPLE_DIR, "wannier90.win"))
    @test mos2_basis.num_wann == 48
    mo_t2g_up = Model.ProjectionGroupConfig(
        "Mo_t2g_up",
        "#1f77b4",
        Int[],
        ["Mo"],
        Int[],
        String[],
        String[],
        ["t2g"],
        "up",
    )
    mo_t2g_spec = Win90Basis.build_projection_spec(mos2_basis, [mo_t2g_up])
    @test length(mo_t2g_spec.groups[1].indices) == 6
    @test mo_t2g_spec.groups[1].atom_count == length(unique(orb.site_label for orb in mos2_basis.orbitals if orb.species == "Mo"))
    s_p = Model.ProjectionGroupConfig(
        "S_p",
        "#2ca02c",
        Int[],
        ["S"],
        Int[],
        String[],
        String[],
        ["p"],
        "any",
    )
    s_p_spec = Win90Basis.build_projection_spec(mos2_basis, [s_p])
    @test length(s_p_spec.groups[1].indices) == 24
    @test s_p_spec.groups[1].atom_count == length(unique(orb.site_label for orb in mos2_basis.orbitals if orb.species == "S"))

    mktempdir() do dir
        spinor_win = joinpath(dir, "spinor_layout.win")
        write(spinor_win, """
num_wann = 8
spinors = .true.

begin unit_cell_cart
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
end unit_cell_cart

begin atoms_frac
X 0.0 0.0 0.0
X 0.5 0.0 0.0
end atoms_frac

begin projections
X : s; pz
end projections
""")
        spinor_qe_basis = Win90Basis.read_win_basis(spinor_win; spin_layout=:qe)
        @test [(orb.site_label, orb.orbital, orb.spin) for orb in spinor_qe_basis.orbitals] == [
            ("X1", "s", :up),
            ("X1", "s", :dn),
            ("X1", "pz", :up),
            ("X1", "pz", :dn),
            ("X2", "s", :up),
            ("X2", "s", :dn),
            ("X2", "pz", :up),
            ("X2", "pz", :dn),
        ]
        spinor_vasp_basis = Win90Basis.read_win_basis(spinor_win; spin_layout=:vasp544)
        @test [(orb.site_label, orb.orbital, orb.spin) for orb in spinor_vasp_basis.orbitals] == [
            ("X1", "s", :up),
            ("X1", "pz", :up),
            ("X2", "s", :up),
            ("X2", "pz", :up),
            ("X1", "s", :dn),
            ("X1", "pz", :dn),
            ("X2", "s", :dn),
            ("X2", "pz", :dn),
        ]
        x_s_down = Model.ProjectionGroupConfig(
            "X_s_down",
            "#1f77b4",
            Int[],
            ["X"],
            Int[],
            String[],
            ["s"],
            String[],
            "dn",
        )
        source_spec = Win90Basis.build_projection_spec(spinor_vasp_basis, [x_s_down])
        @test source_spec.groups[1].indices == [5, 7]

        vasp_win_projection = Model.ProjectionConfig(
            true,
            :win_groups,
            ["X_s_down"],
            :colorbar,
            "viridis",
            String[],
            9.0,
            1.0,
            "weights.dat",
            "projected_bands.png",
            "pdos.dat",
            "pdos.png",
            "combined.png",
            spinor_win,
            [x_s_down],
        )
        vasp_cfg = Model.RunConfig(
            base_cfg.mode,
            base_cfg.files,
            base_cfg.output,
            base_cfg.dos,
            base_cfg.energy,
            base_cfg.plot,
            base_cfg.band_plot,
            base_cfg.dos_plot,
            base_cfg.combined_plot,
            Model.SpinConfig(true, :vasp544, ("#1f77b4", "#d62728")),
            vasp_win_projection,
            base_cfg.hermiticity_tol,
            base_cfg.verbose,
        )
        canonical_spec = Projection.build_projection_spec(vasp_cfg, 8)
        @test canonical_spec.groups[1].indices == [2, 6]
        @test canonical_spec.groups[1].indices != [3, 4]

        win_path = joinpath(dir, "ru_d.win")
        write(win_path, """
num_wann = 5

begin unit_cell_cart
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
end unit_cell_cart

begin atoms_frac
Ru 0.0 0.0 0.0
end atoms_frac

begin projections
Ru : d
end projections
""")
        ru_d_basis = Win90Basis.read_win_basis(win_path)
        ru_t2g_from_d = Model.ProjectionGroupConfig(
            "Ru_t2g",
            "#1f77b4",
            Int[],
            ["Ru"],
            Int[],
            String[],
            ["dxz", "dxy", "dyz"],
            String[],
            "any",
        )
        @test Win90Basis.build_projection_spec(ru_d_basis, [ru_t2g_from_d]).groups[1].indices == [2, 3, 5]
        @test occursin("not defined in win", error_message(() -> Win90Basis.build_projection_spec(ru_d_basis, [ru_p])))

        uneven_win = joinpath(dir, "uneven_ru.win")
        write(uneven_win, """
num_wann = 3

begin unit_cell_cart
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
end unit_cell_cart

begin atoms_frac
Ru 0.0 0.0 0.0
Ru 0.5 0.0 0.0
end atoms_frac

begin projections
f=0.0,0.0,0.0:dxz,dxy
f=0.5,0.0,0.0:dxz
end projections
""")
        uneven_basis = Win90Basis.read_win_basis(uneven_win)
        uneven_group = Model.ProjectionGroupConfig(
            "Ru_dxz_dxy",
            "#1f77b4",
            Int[],
            ["Ru"],
            Int[],
            String[],
            ["dxz", "dxy"],
            String[],
            "any",
        )
        @test_logs (:warn, r"selected site Ru2 lacks requested orbitals") Win90Basis.build_projection_spec(
            uneven_basis,
            [uneven_group],
        )
        uneven_spec = @test_logs (:warn, r"selected site Ru2 lacks requested orbitals") Win90Basis.build_projection_spec(
            uneven_basis,
            [uneven_group],
        )
        @test uneven_spec.groups[1].atom_count == 2

        p_win = joinpath(dir, "p_basis.win")
        write(p_win, """
num_wann = 3

begin unit_cell_cart
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
end unit_cell_cart

begin atoms_frac
A 0.0 0.0 0.0
end atoms_frac

begin projections
A : p
end projections
""")
        local_axes = LocalAxisRotation.parse_axes([
            Any["A1", [0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]],
        ])
        rotation_projection = Model.ProjectionConfig(
            true,
            :win_groups,
            ["A_px"],
            :colorbar,
            "viridis",
            String[],
            9.0,
            1.0,
            "weights.dat",
            "projected_bands.png",
            "pdos.dat",
            "pdos.png",
            "combined.png",
            p_win,
            [
                Model.ProjectionGroupConfig("A_px", "blue", Int[], ["A"], Int[], String[], ["px"], String[], "any"),
            ],
            Model.ProjectionBasisRotationConfig(true, local_axes, false, 1.0e-8),
        )
        rotation_cfg = Model.RunConfig(
            base_cfg.mode,
            base_cfg.files,
            base_cfg.output,
            base_cfg.dos,
            base_cfg.energy,
            base_cfg.plot,
            base_cfg.band_plot,
            base_cfg.dos_plot,
            base_cfg.combined_plot,
            base_cfg.spin,
            rotation_projection,
            base_cfg.hermiticity_tol,
            base_cfg.verbose,
        )
        transform = Projection.build_basis_rotation_transform(rotation_cfg, 3)
        rotated = Projection.projection_weights(
            OrbitalProjection.ProjectionSpec(
                [OrbitalProjection.ProjectionGroup("A_px", [2], "blue")],
                3,
            ),
            Matrix{ComplexF64}(I, 3, 3);
            basis_transform=transform,
        )
        @test rotated[:, 1] ≈ [0.0, 0.0, 1.0]

        input_path = joinpath(dir, "input.rotation.toml")
        hr_path = joinpath(dir, "empty_hr.dat")
        write(hr_path, "")
        write(input_path, """
[run]
mode = "dos"
hr = "empty_hr.dat"

[projection]
enabled = true
mode = "win_groups"
win = "p_basis.win"

[projection.basis_rotation]
enabled = true
local_axes = [["A1", [0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]]]

[[projection.groups]]
label = "A_px"
species = ["A"]
orbitals = ["px"]
""")
        parsed_rotation_cfg = InputIO.read_input(input_path)
        @test parsed_rotation_cfg.projection.basis_rotation.enabled
        @test only(parsed_rotation_cfg.projection.basis_rotation.local_axes).site == "A1"
        @test only(parsed_rotation_cfg.projection.basis_rotation.local_axes).source_x == [0.0, 1.0, 0.0]
    end
end
