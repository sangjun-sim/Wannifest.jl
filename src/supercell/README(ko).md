# supercell

언어: [English](README.md) | 한국어

## 이 모듈이 하는 일

`supercell`은 POSCAR/CONTCAR structure input에서 structural supercell을 만듭니다.

입력 POSCAR/CONTCAR basis에 integer supercell matrix를 그대로 적용해 새
POSCAR를 씁니다. Spglib symmetry summary를 출력하고 symmetry mapping을
검증합니다.

## 언제 사용하나

다음이 필요할 때 `supercell`을 사용합니다.

- structure를 더 큰 POSCAR로 확장,
- supercell construction 전후의 symmetry summary 확인,
- atom-symmetry mapping check를 위한 mismatch diagnostic 작성.

## 엔트리포인트

메인 엔트리:

- `../main.jl` (`spectrum/main.jl`)

모듈 CLI:

- `CLI.jl`

핵심 구현:

- `Service.jl`

## 주요 파일

- `CLI.jl`: argument를 파싱하고 symmetry/check summary를 출력합니다.
- `InputIO.jl`: CLI matrix/path argument와 고정 command 기본값으로
  `RunConfig`를 만듭니다.
- `Service.jl`: input structure를 로드하고 원래 input basis에서 supercell을 만들고 POSCAR output을 쓰며, 선택적 validation을 실행합니다.
- `Core.jl`: supercell submodule을 로드합니다.
- `Symmetry.jl`: Spglib-backed dataset summary helper.
- `Transform.jl`: integer-matrix supercell construction.
- `Validate.jl`: 선택적 symmetry-mapping validation과 mismatch output.

## 동작 방식

상위 흐름:

```text
CLI -> InputIO -> Service -> load POSCAR/CONTCAR structure
                         -> optional Spglib dataset summary
                         -> Transform.build_supercell
                         -> write POSCAR
                         -> Validate checks
```

structural basis는 항상 원래 입력 POSCAR/CONTCAR basis입니다.
CLI는 `use_symmetry=true`, `validate=true`, `symprec=1.0e-5`,
`angle_tolerance=-1.0`, `digits=12`를 고정합니다.

## 입력과 출력

최소 CLI 예:

```bash
julia main.jl supercell 2 2 1
julia main.jl supercell 2 0 0, 0 2 0, 0 0 1
```

명령은 `POSCAR`를 읽고 `POSCAR.supercell`을 쓰며, mismatch
table은 output POSCAR 옆의 `atom_symm_mismatch.dat`로 씁니다.

출력:

- supercell POSCAR,
- stdout symmetry/check summary,
- 선택적 atom-symmetry mismatch table.

## 의존성

내부 의존성:

- `../shared/CrystalCell.jl`
- `../shared/CellConventions.jl`
- `../shared/PoscarIO.jl`
- `Core.jl`
- `Symmetry.jl`
- `Transform.jl`
- `Validate.jl`

외부 의존성:

- `Spglib`
- Julia `LinearAlgebra`

## 최소 코드 예

CLI 사용:

```bash
julia main.jl supercell 2 2 1
```

Service layer 사용:

```julia
include("spectrum/supercell/CLI.jl")

matrix = Main.SupercellCLI.InputIO.parse_matrix_args(["2", "2", "1"])
cfg = Main.SupercellCLI.InputIO.default_config(matrix)
result = Main.SupercellCLI.Service.run(cfg)
println(result.output_path)
```

## 실행과 테스트

모듈 실행:

```bash
julia main.jl supercell 2 2 1
```

CLI smoke-check:

```bash
julia --project=spectrum spectrum/main.jl supercell --help
```

## 흔한 실수

- `supercell`은 `--input`, `--output` 등 추가 argument를 받지 않습니다.
- 명령은 항상 `POSCAR`를 읽습니다. 파일이 없으면 에러가 납니다.
- full matrix row 구분에는 semicolon이 아니라 comma를 씁니다.
- matrix entry 세 개는 diagonal 값으로 해석합니다. full matrix는 세 개의
  comma-separated row이고 각 row는 integer entry 세 개여야 합니다.
- CLI에서는 validation이 항상 켜져 있고 mismatch table을 씁니다.

## 읽는 순서

1. `CLI.jl`
2. `InputIO.jl`
3. `Service.jl`
4. `Symmetry.jl`
5. `Transform.jl`
6. `Validate.jl`
