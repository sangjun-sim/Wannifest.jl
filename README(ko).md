# Wannifest.jl

Wannifest.jl은 Wannier tight-binding 모델과 결정 구조 입력을 다루기 위한
Julia toolkit입니다. 하나의 command-line entry point에서 band 분석,
k-plane energy map, flux term 생성, structural supercell 생성, supercell
Hamiltonian 생성을 실행할 수 있습니다.

이 package는 Wannier90 스타일 `hr.dat`, `KPOINTS`, POSCAR/CONTCAR 구조 파일,
`wannier90.win` projection 정보, TOML 입력 파일을 사용하는 workflow를
대상으로 합니다. 명시적인 입력 schema, 재현 가능한 출력 데이터, 독립적으로
실행 가능한 command 단위 도구를 지향합니다.

## 기능

- `band`: band structure, DOS, projection, OAM, SAM workflow.
- `contour`: 2D k-plane energy surface와 contour/surface plot 생성.
- `flux`: Wannier Hamiltonian에 directed complex flux term 추가, optional 3D
  HTML visualization 및 diagnostic 출력.
- `supercell`: POSCAR/CONTCAR structural supercell 생성, symmetry summary 및
  validation 출력.
- `superham`: primitive Hamiltonian을 supercell Hamiltonian으로 확장하고,
  optional orbital center 및 Wigner-Seitz vector 정보를 처리.

## 설치

Repository root에서 Julia environment를 instantiate합니다.

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## 사용법

Command 목록을 확인합니다.

```bash
julia --project=. main.jl --help
```

Plot 없이 band 계산을 실행합니다.

```bash
julia --project=. main.jl band --input examples/band/rucl3/input.toml --no-plot
```

k-plane energy surface를 sampling합니다.

```bash
julia --project=. main.jl contour --input examples/contour/graphene/input.toml --no-plot
```

현재 directory의 `POSCAR`에서 structural supercell을 생성합니다.

```bash
julia --project=. main.jl supercell 2 0 0, 0 2 0, 0 0 1
```

Directed flux term을 추가하고 Hamiltonian 출력만 생성합니다.

```bash
julia --project=. main.jl flux --input examples/flux/graphene/input.toml --no-html --no-diagnostic
```

TOML 입력에서 supercell Hamiltonian을 생성합니다.

```bash
julia --project=. main.jl superham --input path/to/input.superham.toml
```

TOML 기반 command는 input path 또는 top-level TOML namespace에서 command를
추론할 수도 있습니다.

```bash
julia --project=. main.jl --input examples/band/rucl3/input.toml
julia --project=. main.jl examples/contour/graphene/input.toml
```

## 입력과 출력

대부분의 command는 `[band.run]`, `[contour.plane]`, `[flux.terms]`,
`[superham.supercell]` 같은 namespaced TOML table을 사용합니다. 경로는
command-line override가 없는 한 input file 위치를 기준으로 해석됩니다.

생성된 data와 plot file은 기본적으로 command별 `outputs/` directory 아래에
기록됩니다. 예제 출력은 git에서 무시되므로 예제를 다시 실행해도 working tree가
불필요하게 더러워지지 않습니다.

## Repository 구조

- `main.jl`: executable command-line entry point.
- `src/Wannifest.jl`: command dispatch 및 TOML command inference.
- `src/shared/`: Wannier, lattice, spin-layout, projection, parsing 공통 utility.
- `src/band/`, `src/contour/`, `src/flux/`, `src/supercell/`,
  `src/superham/`: command implementation.
- `examples/`: 실행 가능한 example input과 reference data.
- `test/`: integration 및 module test.

## 테스트

```bash
julia --project=. --startup-file=no test/runtests.jl
```
