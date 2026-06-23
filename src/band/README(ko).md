# band

언어: [English](README.md) | 한국어

## 이 모듈이 하는 일

`band`는 Wannier 형식의 `hr.dat`에서 다음을 계산합니다.

- band structure,
- density of states, DOS,
- orbital-projected band weight와 projected DOS, PDOS,
- band path의 projected local atomic OAM,
- band path의 spin angular momentum expectation,
- 명시적 spin layout을 요청한 경우 spin-colored plot과 spin-resolved DOS,
- 위 출력들의 지원되는 조합.

명시적 opt-in `wsvec` correction, POSCAR/CONTCAR 또는 `wannier90.win`
lattice input, `Plots` lazy plotting, strict TOML validation을 지원합니다.

## 언제 사용하나

다음이 필요할 때 `band`를 사용합니다.

- `KPOINTS` 파일에서 band path plot 만들기,
- uniform mesh에서 DOS curve 계산하기,
- band path나 DOS에서 orbital character 확인하기,
- 서로 다른 Hamiltonian의 band 비교하기,
- 외부 plotting을 위한 data file 생성하기.

## 엔트리포인트

메인 엔트리:

- `../main.jl`

모듈 CLI:

- `CLI.jl`

핵심 구현:

- `Service.jl`

## 주요 파일

- `CLI.jl`: `--input`, `--no-plot`을 파싱합니다.
- `InputIO.jl`: TOML을 읽고 run, output, DOS, plotting, spin, projection,
  OAM option을 검증합니다.
- `Service.jl`: band와 DOS 실행의 주요 orchestration.
- `../shared/LatticeIO.jl`: POSCAR/CONTCAR 또는 `wannier90.win` lattice 입력을 읽습니다.
- `KPath.jl`: `KPOINTS`를 파싱하고 보간된 path를 생성합니다.
- `Execution.jl`: scoped execution 동안 BLAS thread 수를 임시로 설정합니다.
- `Bands.jl`: band-structure post-processing.
- `Dos.jl`: DOS mesh generation과 broadening.
- `projection/Model.jl`: projection config와 result type.
- `projection/Input.jl`: `[band.projection]` parser와 기본 output path.
- `projection/InputGroups.jl`: compact/table projection group parser.
- `PlotInput.jl`: `[band.plot]` parser와 명시적 plot target 검증.
- `projection/Core.jl`: projection spec과 projection output을 만듭니다.
- `observables/Oam.jl`: projection basis metadata에서 projected atomic OAM
  context를 만듭니다.
- `observables/Sam.jl`: spinful basis metadata에서 spin angular momentum
  context를 만듭니다.
- `observables/Output.jl`: long-format OAM과 SAM data를 씁니다.
- `observables/Plot.jl`: OAM/SAM color band structure를 그립니다.
- `../shared/SpinLayout.jl`: 공유 명시적 spin-layout index, spin coloring,
  spin summary text.
- `Output.jl`: band, DOS, projection weight, PDOS data writer.
- `PlotService.jl`: plotting lazy loader와 dispatcher.
- `PlotBands.jl`: `Plots` plotting backend. Plotting이 요청될 때만 load됩니다.
- `../shared/OrbitalProjection.jl`: group 정의와 eigenvector weight 계산.
- `../shared/AtomicOam.jl`: local atomic `Lx`, `Ly`, `Lz`, `L2` operator.
- `../shared/Win90Basis.jl`: `wannier90.win` projection block을 읽어 selector 기반 group을 만듭니다.
- `../shared/WannierEigensystem.jl`: `H(k)`를 diagonalize합니다.

## 동작 방식

상위 흐름:

```text
CLI -> InputIO -> Service -> layout-aware WannierHrIO.read_hr
                         -> optional layout-aware WannierWsvecIO.read_wsvec + validation
                         -> optional LatticeIO.read_lattice
                         -> KPath / Dos mesh generation
                         -> WannierEigensystem.solve_kpoint
                         -> explicit spin layout split when enabled
                         -> optional projection weights / PDOS
                         -> optional projected atomic OAM
                         -> Output data writers
                         -> optional PlotService plots
```

일반적인 compute path는 세 가지입니다.

- band path: `KPOINTS`를 파싱하고 path를 보간한 뒤 eigenvalue를 풀어 band data를 씁니다.
- DOS path: regular mesh를 만들고 많은 `k` point를 푼 뒤 eigenvalue를 broaden해서 DOS data를 씁니다.
- projection path: 같은 diagonalization에서 나온 sorted eigenvector를 사용해 group index의 `abs2` 합을 계산하고 별도 projection output을 씁니다.
- OAM path: 같은 diagonalization의 sorted eigenvector로 projected local atomic
  OAM을 계산하고 별도 long-format data를 씁니다.

## 입력과 출력

최소 band 예:

```toml
[band.run]
mode = "bands"
# 생략 시 기본값:
# hr = "wannier90_hr.dat"
# wsvec 생략 또는 wsvec = "" 는 wsvec correction을 끕니다.
# wsvec = "wannier90_wsvec.dat" 처럼 명시해야 사용합니다.
# kpoints = "KPOINTS"
# structure = "POSCAR"
# hermiticity_tol = 1.0e-8
# verbose = true

[band.energy]
shift = 0.0
```

최소 DOS 예:

```toml
[band.run]
mode = "dos"
hr = "wannier90_hr.dat"
structure = "POSCAR"

[band.dos]
mesh = "50x50x50"
shift = [0.0, 0.0, 0.0]
sigma = 0.05
npts = 2001
emin = -10
emax = 10

[band.energy]
shift = 0.0
```

Combined band/DOS label 예:

```toml
[band.combined_plot]
dos_width_ratio = 0.30
dos_xlabel = "DOS (states/eV)"  # 회전된 DOS panel의 x label
dos_ylabel = ""                 # 회전된 DOS panel의 y label
```

Plot 파일은 `[band.plot].targets`로 명시한 경우에만 씁니다. `--no-plot`은
targets가 있어도 모든 plotting을 비활성화합니다.

```toml
[band.plot]
targets = ["band", "dos", "combined", "fatband", "pdos", "fatband_pdos", "oam", "sam"]
size = [900, 600]
energy_range = [-3.0, 3.0]
font_size = 18
```

허용 target은 `band`, `dos`, `combined`, `fatband`, `pdos`, `fatband_pdos`,
`oam`, `sam`입니다. `targets`를 생략하거나 빈 배열로 두면 data 파일만 쓰고 plot
파일은 쓰지 않습니다. Projection target은 `[band.projection].enabled = true`,
`oam` target은 `[band.oam].enabled = true`, `sam` target은
`[band.sam].enabled = true`가 필요합니다.

DOS plot에는 center-of-mass marker를 dotted line으로 함께 그립니다. Standalone
DOS plot에서는 energy 축의 vertical line이고, combined plot에서는 회전된 DOS
panel의 horizontal line입니다. Spin-resolved DOS는 `dos_up + dos_down` total
center를 사용하며, PDOS plot은 group별 center marker를 해당 group color로
표시합니다.

Spin-colored band와 spin-resolved DOS 예:

```toml
[band.spin]
enabled = true
layout = "qe"  # "qe" 또는 "vasp544"
colors = ["#1f77b4", "#d62728"]
```

`layout = "qe"`는 QE/interleaved spin ordering(`up, down, up, down, ...`)을
사용합니다. `layout = "vasp544"`는 VASP 5.4.4 block spin ordering
(`all up, then all down`)을 사용합니다. layout은 명시 입력이며, `band`는
Hamiltonian에서 layout을 추론하지 않습니다. VASP block order의 HR, wsvec,
projection index는 diagonalization 전에 내부 QE/interleaved order로
canonicalize됩니다.
Band plot은 모든 layout에서 정렬된 eigenvalue column 순서를 그대로 따릅니다.
Spin coloring과 spin-resolved DOS는 입력 canonicalization 이후의 내부
canonical order(`up, down, up, down, ...`)를 사용합니다.

명시적 Wannier index 기반 orbital projection 예:

```toml
[band.projection]
enabled = true
mode = "index_groups"
plot_style = "colorbar"                 # "colorbar" 또는 "empty_circle"
colorbar_colormap = "viridis"           # Plots.jl color gradient 이름
groups = [
  ["Ru_t2g", [1, 2, 3], "#1f77b4"],
  ["Cl_p", [4, 5, 6], "#2ca02c"],
]
```

`wannier90.win` 기반 orbital projection 예:

```toml
[band.projection]
enabled = true
mode = "win_groups"
win = "wannier90.win"
groups = [
  ["Mo_sp", ["Mo"], ["s", "p"], "red"],
  ["Mo_t2g", ["Mo"], ["t2g"], "blue"],
]
```

지원하는 projection mode:

- `index_groups`: group을 1-based Wannier orbital index로 직접 지정합니다.
- `win_groups`: `wannier90.win`의 `begin projections` block에서 selector로 group을 만듭니다.
- `color_group`은 선택 사항입니다. 생략하면 projected band plot은 `groups`에
  적힌 모든 group을 사용합니다. 일부 group만 그리고 싶을 때만
  `color_group = ["label_a"]`처럼 지정합니다.
- `groups`는 compact row 형식도 받을 수 있습니다. `index_groups`에서는
  `[label, indices, color]`, `win_groups`에서는 `[label, species, orbitals, color]`를 씁니다.
  `site_labels`, 숫자 `sites`, `spin` 같은 세부 selector가 필요하면
  `[[band.projection.groups]]` table 형식도 계속 사용할 수 있습니다.

지원하는 projected-band plot style:

- `plot.targets`로 요청한 경우, projected band plot은 선택된 projection
  group을 overlay하고 PDOS plot은 항상 전체 group을 표시합니다.
- `plot_style = "colorbar"`: 고정 크기 marker를 선택된 각 group weight로 색칠합니다. 여러 group을 선택하면 group별 색 gradient를 사용합니다.
- `plot_style = "empty_circle"`: 빈 원을 그리고, 선택된 각 group weight에 비례해 원 크기를 정합니다. Stroke color는 group color를 사용합니다. `circle_max_size`, `circle_stroke_width`로 marker geometry를 조정할 수 있습니다.
- `colorbar_colormap = "plasma"`로 colorbar mode의 named gradient를 바꿀 수 있습니다. 직접 gradient를 지정하려면 `colorbar_colors = ["white", "#1f77b4"]`를 씁니다.

Projection data output은 기존 band/DOS output과 별도입니다. 기존
`bands.dat`, `dos.dat`는 바꾸지 않습니다.
PDOS curve는 projection source에 site metadata가 있는 경우(`win_groups`) group별 atom 수로 나눈 per-atom 값입니다. `index_groups`는 atom metadata가 없으므로 atom count를 1로 둡니다.

Projected atomic OAM은 같은 projection basis metadata를 사용합니다. Band path
run에서만 동작하며, `[band.projection]`이 `mode = "win_groups"`로 metadata를 제공해야 합니다.

```toml
[band.oam]
enabled = true
orbitals = [
  ["Ru1", "t2g"],
  ["Cl1", "p"],
]
degeneracy_tol = 1.0e-4
plot_components = ["Lz"]
```

각 `orbitals` row는 `[site_label, shell]`입니다. 첫 값은 site label로만
해석하며 species fallback은 쓰지 않습니다. 첫 버전은 shell token
`s`, `p`, `d`, `t2g`, `eg`를 지원합니다. `plot_components` 기본값은
`["Lz"]`이며 `Lx`, `Ly`, `Lz`, `L_norm`, `L2`를 지정할 수 있습니다.
OAM은 `[band.output] bands_data`에서 파생한 `bands_oam.dat`를 씁니다.
OAM color band plot은 `plot.targets`에 `oam`이 있을 때만 쓰며,
`[band.output] bands_plot`에서 파생한 `bands_oam_lz.png` 같은 경로를
사용합니다. 이 파생 파일들의 TOML 출력 경로는 받지 않습니다. 이 값은
선택한 Wannier basis의 projected local atomic OAM이며 Berry-phase orbital
magnetization이 아닙니다. 정확한 degeneracy에서는 per-band OAM이
gauge-dependent할 수 있고, 첫 버전은 경고만 내며 degenerate-cluster trace는
쓰지 않습니다.

Spin angular momentum expectation은 같은 projection basis metadata를 사용하며
`[band.sam]`으로 켭니다. 명시적 `[band.spin] layout`이 필요하지만
`spin.enabled = true`는 필수가 아닙니다.

```toml
[band.spin]
layout = "vasp544" # 또는 "qe"

[band.sam]
enabled = true
degeneracy_tol = 1.0e-4
plot_components = ["Sz"]
```

SAM은 `[band.output] bands_data`에서 파생한 `bands_sam.dat`를 쓰며 `Sx`,
`Sy`, `Sz`, `S_norm`, `S2`를 기록합니다. `plot_components` 기본값은
`["Sz"]`이며 `Sx`, `Sy`, `Sz`, `S_norm`, `S2`를 지정할 수 있습니다. SAM color
band plot은 `plot.targets`에 `sam`이 있을 때만 쓰며 `[band.output] bands_plot`에서
파생한 `bands_sam_sz.png` 같은 경로를 사용합니다. Pauli block은 complete
`(site, orbital, up/dn)` pair가 필요하므로 nonmagnetic 또는 unpolarized basis
metadata는 error입니다.

기본 data output과 plot destination은 입력 TOML 디렉터리 기준으로
해석됩니다. Plot destination은 대응하는 `plot.targets` 항목이 있을 때만
사용됩니다.

| Table/key | 기본값 |
| --- | --- |
| `[band.output] bands_data` | `outputs/data/bands.dat` |
| `[band.output] dos_data` | `outputs/data/dos.dat` |
| `[band.output] bands_plot` | `outputs/plots/bands.png` |
| `[band.output] dos_plot` | `outputs/plots/dos.png` |
| `[band.output] combined_plot` | `outputs/plots/band_dos.png` |
| `[band.projection] weights_data` | `outputs/data/projection/band_projection_weights.dat` |
| `[band.projection] pdos_data` | `outputs/data/projection/pdos.dat` |
| `[band.projection] projected_bands_plot` | `outputs/plots/projection/bands_projected.png` |
| `[band.projection] pdos_plot` | `outputs/plots/projection/pdos.png` |
| `[band.projection] projected_combined_plot` | `outputs/plots/projection/band_pdos.png` |
| `[band.oam] 파생 data` | `outputs/data/bands_oam.dat` |
| `[band.oam] 파생 plot` | `outputs/plots/bands_oam_lz.png` |
| `[band.sam] 파생 data` | `outputs/data/bands_sam.dat` |
| `[band.sam] 파생 plot` | `outputs/plots/bands_sam_sz.png` |

실행 예:

```bash
julia main.jl band --input input.band.toml
```

```bash
julia main.jl band --input input.band.toml --no-plot
```

SAM 예:

```bash
julia spectrum/main.jl spectrum/examples/sam/bi2se3/input.toml --no-plot
```

## 의존성

내부 의존성:

- `../shared/CrystalCell.jl`
- `../shared/PoscarIO.jl`
- `../shared/LatticeIO.jl`
- `../shared/WannierHrIO.jl`
- `../shared/WannierWsvecIO.jl`
- `../shared/WannierWsvecGenerate.jl`
- `../shared/WannierKspace.jl`
- `../shared/WannierEigensystem.jl`
- `../shared/OrbitalProjection.jl`
- `../shared/AtomicOam.jl`
- `../shared/AtomicSpin.jl`
- `../shared/Win90Basis.jl`
- `projection/Model.jl`
- `projection/InputGroups.jl`
- `projection/Input.jl`
- `PlotInput.jl`
- `projection/Core.jl`
- `InputIO.jl`
- `KPath.jl`
- `Execution.jl`
- `observables/Oam.jl`
- `observables/Sam.jl`
- `observables/Output.jl`
- `observables/Plot.jl`
- `../shared/SpinLayout.jl`
- `Bands.jl`
- `Dos.jl`
- `Output.jl`
- `PlotService.jl`
- `PlotBands.jl`

외부 의존성:

- `Plots`
- Julia `Printf`
- Julia threading support

## 최소 코드 예

CLI 사용:

```bash
julia main.jl band --input input.band.toml --no-plot
```

Service layer 사용:

```julia
include("spectrum/band/CLI.jl")

cfg = Main.BandCLI.InputIO.read_input("spectrum/examples/graphene/input.toml")
result = Main.BandCLI.Service.run(cfg; make_plot=false)
println(result.hermiticity_ok)
```

## 실행과 테스트

모듈 실행:

```bash
julia main.jl band --input input.band.toml --no-plot
```

테스트 실행:

```bash
julia --project=spectrum --startup-file=no spectrum/test/runtests.jl
```

## 흔한 실수

- `wsvec` correction은 opt-in입니다. `run.wsvec` 생략과 `wsvec = ""`는 둘 다 비활성화이며, 비어 있지 않은 경로는 실제 파일이어야 하고 모든 `(R, i, j)` entry를 포함해야 합니다.
- Cartesian `KPOINTS`는 좌표를 input cell의 reciprocal basis로 reduced coordinate화하기 위해 `run.structure`가 필요합니다. 이 값은 POSCAR/CONTCAR 또는 `wannier90.win`일 수 있습니다.
- Reduced `KPOINTS`는 `hr.dat`가 쓰는 input cell의 reciprocal basis에서 직접 해석됩니다. HR `R` vector는 자동으로 재기저변환하지 않습니다. `run.structure`가 `hr.dat`와 같은 cell convention이 아니면 K-path 거리축이 왜곡될 수 있습니다.
- lattice 없이 reduced-coordinate path distance를 의도적으로 쓰려는 경우에만
  `run.structure = ""`를 설정합니다.
- Projection weight는 raw eigenvector에서 계산합니다. `energy.shift`는 eigenvalue와 DOS energy에만 적용됩니다.
- `win_groups`에서는 사용자가 요청한 orbital을 `wannier90.win`과 대조합니다. 예를 들어 win에는 `Ru:d`만 있는데 `Ru:p`를 요청하면 빈 group을 만들지 않고 에러를 냅니다.
- 정확한 degeneracy에서는 개별 band의 projection weight가 eigenvector gauge에 의존할 수 있습니다. 단일 band 색보다 degenerate subspace의 group 합이 더 안정적입니다.
- 정확한 degeneracy에서는 개별 band OAM도 eigenvector gauge에 의존할 수
  있습니다. `band.oam`은 per-band 값과 경고를 내며 cluster trace는 쓰지 않습니다.
- `band.oam`은 Wannier basis의 projected local atomic observable이며 orbital
  magnetization 계산이 아닙니다.
- 알 수 없는 `[band.*]` table과 지원하지 않는 key는 early failure가 납니다.
  실험적 note는 `[band]` namespace 밖에 둡니다.
- Plot generation은 `Plots`에 의존하므로, 계산은 성공해도 환경 문제 때문에 plotting만 실패할 수 있습니다.
- `energy.shift = 0.0`은 raw eigenvalue를 의미하며, Fermi level로 자동 shift하지 않습니다.

## 읽는 순서

1. `CLI.jl`
2. `InputIO.jl`
3. `PlotInput.jl`
4. `Service.jl`
5. `KPath.jl`
6. `Bands.jl`
7. `Dos.jl`
8. `projection/Core.jl`
9. `Output.jl`
10. `PlotService.jl`
11. `../shared/SpinLayout.jl`
12. `PlotBands.jl`
13. `observables/Plot.jl`
14. `../shared/WannierEigensystem.jl`
15. `../shared/OrbitalProjection.jl`
16. `../shared/AtomicOam.jl`
17. `../shared/AtomicSpin.jl`
18. `observables/Sam.jl`
19. `../shared/Win90Basis.jl`
