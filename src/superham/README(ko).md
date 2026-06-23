# superham

언어: [English](README.md) | 한국어

## 이 모듈이 하는 일

`superham`은 primitive-cell `hr.dat`에서 supercell Hamiltonian을 만듭니다.
구조 lattice는 입력 파일 그대로 사용합니다. `hr.dat`의 `R` vector와 manual center는 같은 input-cell basis 기준이어야 합니다.

다음과 같은 geometry information과 Wigner-Seitz correction도 붙일 수 있습니다.

- `wsvec` table,
- `centres.xyz`의 Wannier center,
- 사용자가 제공한 manual center,
- atom-based surrogate center.

CLI는 요청한 `k` point에서 folded spectrum을 check하고 normalized supercell
`hr.dat`를 쓸 수 있습니다.

## 언제 사용하나

다음이 필요할 때 `superham`을 사용합니다.

- primitive Hamiltonian을 더 큰 supercell Hamiltonian으로 확장,
- geometry-aware workflow를 위해 center information 유지 또는 부착,
- tight-binding model이 supercell folding에서 어떻게 바뀌는지 테스트,
- supercell용 새 `hr.dat` 작성.

## 엔트리포인트

메인 엔트리:

- `../main.jl`

모듈 CLI:

- `CLI.jl`

핵심 구현:

- `Service.jl`

## 주요 파일

- `CLI.jl`: command-line interface.
- `InputIO.jl`: TOML에서 file, geometry, supercell option을 읽습니다.
- `Service.jl`: 주요 orchestration.
- `HrIO.jl`: 입력 `hr.dat`를 읽습니다.
- `CenterIO.jl`: center table을 로드하거나 만듭니다.
- `SupercellHam.jl`: 확장된 Hamiltonian을 구성합니다.
- `Kspace.jl`: `k`에서 eigenvalue를 평가합니다.
- `Validate.jl`: Hermiticity와 folding check.
- `WsvecIO.jl`: 선택적 Wigner-Seitz correction 지원과 generated `wsvec`
  entrypoint.
- `../shared/HrFormat.jl`: normalized supercell `hr.dat` output writer.
- `../test/runtests.jl`: superham helper를 포함하는 shared spectrum test.

## 동작 방식

상위 흐름:

```text
CLI -> InputIO -> Service -> PoscarIO.read_poscar
                         -> HrIO.read_hr
                         -> optional WsvecIO.attach_wsvec
                         -> optional CenterIO.centers_from_config
                         -> SupercellHam.build_supercell_model
                         -> Validate hermiticity / size / folded spectrum
                         -> optional normalized hr.dat writer
```

이 모듈은 여러 geometry mode를 지원합니다.

- `none`: orbital-center table을 붙이지 않습니다.
- `wsvec`: wsvec-only geometry policy를 나타냅니다. 실제 `wsvec` attachment는
  비어 있지 않은 `files.wsvec` path가 제어합니다.
- `manual_centers`: `[[superham.geometry.orbitals]]` entry를 fractional
  Wannier center로 사용합니다.
- `atomic_assumption`: 각 orbital entry의 `atom_index`로 structure file의 atom
  position을 approximate center로 복사합니다.

Geometry mode는 orbital-center information을 생략할지, manual로 줄지,
atom에서 근사할지, explicit `wsvec` data와 함께 쓸지를 결정합니다.

## 입력과 출력

최소 입력 예:

```toml
[superham.files]
hr = "wannier90_hr.dat"
structure = "POSCAR"
win = ""
wsvec = ""
centres = ""
output_hr = "wannier90_sc_hr.dat"

[superham.geometry]
mode = "none"
strict = false

[superham.supercell]
matrix = [2, 0, 0, 0, 2, 0, 0, 0, 1]
```

Supercell matrix는 nested 3x3 array로도 쓸 수 있습니다.

```toml
[superham.supercell]
matrix = [
  [2, 0, 0],
  [0, 2, 0],
  [0, 0, 1],
]
```

Manual-centers 예:

```toml
[superham.geometry]
mode = "manual_centers"
manual_num_wann = 2

[[superham.geometry.orbitals]]
label = "A:pz"
center_frac = [0.0, 0.0, 0.0]

[[superham.geometry.orbitals]]
label = "B:pz"
center_frac = [0.3333333333, 0.3333333333, 0.0]
```

Atomic-assumption 예:

```toml
[superham.geometry]
mode = "atomic_assumption"
manual_num_wann = 2

[[superham.geometry.orbitals]]
label = "A:pz"
center_frac = [0.0, 0.0, 0.0]
atom_index = 1

[[superham.geometry.orbitals]]
label = "B:pz"
center_frac = [0.0, 0.0, 0.0]
atom_index = 2
```

실행 예:

```bash
julia main.jl superham --input input.superham.toml
julia main.jl superham --input input.superham.toml --kpoint 0 0 0 --output-hr wannier90_sc_hr.dat
```

출력:

- supercell `hr.dat`,
- geometry source, primitive/supercell orbital count, folded-spectrum
  validation, 선택한 `k` point의 eigenvalue가 포함된 stdout summary.

`[superham.files]`의 `output_hr`는 CLI `--output-hr` override가 없을 때
사용됩니다. 둘 다 없으면 HR file은 쓰지 않고 validation과 summary만
실행합니다.

`files.win`은 configuration object에서 받지만, 현재 service는 lattice와 atom
position을 `files.structure`에서 구성합니다.

## 의존성

내부 의존성:

- `../shared/CrystalCell.jl`
- `../shared/HrFormat.jl`
- `../shared/WannierHrIO.jl`
- `../shared/WannierKspace.jl`
- `../shared/WannierWsvecIO.jl`
- `../shared/WannierWsvecGenerate.jl`
- `../shared/PoscarIO.jl`
- `InputIO.jl`
- `CenterIO.jl`
- `SupercellHam.jl`
- `Validate.jl`

외부 의존성:

- Julia `TOML`
- Julia `LinearAlgebra`

## 최소 코드 예

CLI 사용:

```bash
julia main.jl superham --input input.superham.toml
```

Service layer 사용:

```julia
include("spectrum/superham/CLI.jl")

cfg = Main.SuperhamCLI.InputIO.read_input("input.superham.toml")
result = Main.SuperhamCLI.Service.run(cfg)
println(result.output_hr)
```

## 실행과 테스트

모듈 실행:

```bash
julia main.jl superham --input input.superham.toml
```

테스트 실행:

```bash
julia --project=spectrum --startup-file=no spectrum/test/runtests.jl
```

## 흔한 실수

- structure, HR `R` vector, manual center 좌표는 같은 input-cell basis를 써야 합니다.
- `geometry.mode = "manual_centers"`와 `geometry.mode = "atomic_assumption"`은 `geometry.orbitals`가 필요합니다.
- `manual_num_wann`은 `hr.dat` basis size와 일치해야 합니다.
- `atomic_assumption`은 Wannier center가 실제로 atom-centered일 때만 정확합니다.
  이 mode를 쓰면 코드가 warning을 출력합니다.
- path가 설정되어 있고 파일이 존재하면 `centres.xyz`가 manual 또는 atomic
  center generation보다 우선합니다.
- `files.wsvec = ""`는 `wsvec`를 끕니다. 비어 있지 않은 path는 attach 전에
  coverage validation을 통과해야 합니다.
- `geometry.strict = true`는 geometry-sensitive supercell construction에서
  missing center information에 더 보수적으로 대응합니다.

## 읽는 순서

1. `CLI.jl`
2. `InputIO.jl`
3. `Service.jl`
4. `HrIO.jl`
5. `CenterIO.jl`
6. `SupercellHam.jl`
7. `Validate.jl`
8. `../test/runtests.jl`
