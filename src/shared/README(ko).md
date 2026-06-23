# shared

언어: [English](README.md) | 한국어

## 이 모듈이 하는 일

`shared`는 `spectrum/` 내부 모듈들이 사용하는 building block을 담고 있습니다.

사용자-facing CLI 모듈은 아닙니다. 대신 다음을 제공합니다.

- crystal-cell과 lattice representation,
- POSCAR/CONTCAR와 `wannier90.win` lattice 입력,
- path resolution을 포함한 TOML parsing helper,
- Wannier `hr.dat`와 `wsvec` data type 및 IO,
- HR 계열 block dictionary용 Hermiticity completion 및 pair-symmetry check,
- `H(k)` construction과 dense eigensystem solve,
- Wigner-Seitz vector generation과 validation,
- structural supercell에서 쓰는 Spglib-backed cell convention helper,
- explicit index, sidecar basis file, `wannier90.win` 기반 projection group,
- Wannier-basis observable용 projected local atomic OAM operator,
- OAM과 SOC가 공유하는 atomic angular momentum matrix convention,
- spinful Wannier-basis observable용 spin angular momentum operator,
- spinless-to-spinful HR 및 basis expansion helper,
- supercell export에서 쓰는 HR writing helper.

## 언제 사용하나

다음을 이해하고 싶을 때 `shared`를 읽습니다.

- 모듈들이 `CrystalCell`을 공유하는 방식,
- POSCAR/CONTCAR와 `wannier90.win` lattice input이 표현되는 방식,
- `hr.dat`, normalized hopping, `wsvec` table이 표현되는 방식,
- band, superham, flake가 TOML parsing utility 중복을 피하는 방식,
- projection selector가 Wannier orbital index로 매핑되는 방식,
- orbital label에서 projected atomic OAM matrix를 조립하는 방식,
- complete `up/dn` basis pair에서 spin operator를 조립하는 방식,
- 생성된 `wsvec` data를 사용 전에 검증하는 방식,
- structural supercell generation이 POSCAR/CONTCAR 입력을 standardize하는 방식.

## 엔트리포인트

독립 CLI는 없습니다.

상위 모듈을 통해 helper를 로드합니다.

```julia
include("spectrum/band/CLI.jl")
include("spectrum/superham/Core.jl")
```

각 module loader는 자신에게 필요한 shared file만 include합니다.

## 주요 파일

- `CrystalCell.jl`: 공유 crystal-cell representation과 lattice helper.
- `PoscarIO.jl`: `CrystalCell` 기반 공유 POSCAR reader.
- `LatticeIO.jl`: POSCAR/CONTCAR와 `wannier90.win`의 `unit_cell_cart` block을
  읽는 lattice reader.
- `InputParsing.jl`: 재사용 가능한 TOML helper 함수.
- `CellConventions.jl`: cell convention, lattice transform, Spglib-backed
  standardized-cell helper.
- `WannierTypes.jl`: HR과 WS-vector data를 위한 공유 type.
- `SpinLayout.jl`: 공유 `vasp544`/`qe` spin-layout parsing과 index
  canonicalization.
- `WannierHrIO.jl`: 공유 `hr.dat` reading helper.
- `Hermiticity.jl`: block dictionary의 pair-Hermiticity error 계산과 assertion.
- `HrHermiticity.jl`: Hermitian partner completion과 onsite block cleanup.
- `PairChecks.jl`: 상위 모듈이 쓰는 pair-dictionary validation helper.
- `WannierKspace.jl`: plain 및 `wsvec`-aware `k`-space Hamiltonian evaluation.
- `WannierEigensystem.jl`: Wannier Hamiltonian용 dense eigenvalue/eigenvector solve.
- `WannierWsvecIO.jl`: 공유 `wsvec` handling.
- `WannierWsvecGenerate.jl`: generated `wsvec` construction, validation, writer utility.
- `OrbitalProjection.jl`: projection group data model과 eigenvector weight.
- `AtomicAngularMomentum.jl`: canonical p/d/t2g angular momentum matrix.
- `AtomicOam.jl`: local atomic `Lx`, `Ly`, `Lz`, projected `L2` matrix.
- `AtomicSpin.jl`: spinful basis용 spin `Sx`, `Sy`, `Sz`, `S2` matrix.
- `SpinExpand.jl`: spinless HR/basis expansion과 collinear up/down merge helper.
- `BasisLabelNormalize.jl`: 공유 orbital 및 spin label canonicalization.
- `Win90Basis.jl`: `wannier90.win` projection-basis reader와 selector expansion.
- `HrFormat.jl`: 여러 모듈이 공유하는 HR writing helper.

## 동작 방식

`shared`는 단독으로 실행되지 않습니다.

대신 다른 모듈이 그 일부를 가져다 씁니다.

```text
band     -> LatticeIO, Hermiticity, PairChecks, WannierHrIO, WannierKspace, WannierEigensystem,
            SpinLayout, WannierWsvecIO, WannierWsvecGenerate, OrbitalProjection,
            AtomicOam, Win90Basis
soc      -> HrFormat, HrHermiticity, WannierHrIO, SpinLayout, SpinExpand,
            AtomicAngularMomentum, Win90Basis
superham -> HrFormat, CrystalCell, PoscarIO, SpinLayout, Hermiticity, PairChecks,
            WannierHrIO, WannierKspace, WannierEigensystem, WannierWsvecIO,
            WannierWsvecGenerate
supercell -> CrystalCell, CellConventions, PoscarIO, InputParsing
flake    -> superham/Core.jl을 통해 shared helper 사용
```

`InputParsing.namespaced_root(cfg, "band")`는 `[band.run]` 같은 namespaced
입력과 `[run]` 같은 module-local 입력을 모두 받을 수 있게 합니다.

## 입력과 출력

`shared`에 대한 단일 입력 파일은 없습니다.

lattice reader는 입력 cell basis를 그대로 보존합니다.

```julia
include("spectrum/shared/CrystalCell.jl")
include("spectrum/shared/PoscarIO.jl")
include("spectrum/shared/LatticeIO.jl")

lattice = Main.LatticeIO.read_lattice("POSCAR")
```

`read_lattice`는 POSCAR/CONTCAR 계열 파일과 `begin unit_cell_cart` block을
가진 `.win` 파일을 받습니다. TOML 파일은 lattice input으로 거부합니다.

module input parser가 사용하는 path resolver는 상대 경로를 입력 TOML
디렉터리 기준으로 유지합니다.

```julia
include("spectrum/shared/InputParsing.jl")

path = Main.InputParsing.resolve_path(dirname(abspath("input.toml")), "KPOINTS")
```

Projection helper는 group label, 1-based index, duplicate coverage, group set이
전체 Wannier orbital을 덮는지 여부를 검증합니다.

`WannierWsvecGenerate.assert_wsvec_usable(hr, wsvec)`는 workflow가 `wsvec`를
사용하기 전에 missing, out-of-bounds, malformed, duplicate, one-sided entry를
검증합니다.

## 의존성

내부 의존성:

- `CrystalCell.jl`
- `PoscarIO.jl`
- `LatticeIO.jl`
- `InputParsing.jl`
- `CellConventions.jl`
- `WannierTypes.jl`
- `SpinLayout.jl`
- `WannierHrIO.jl`
- `WannierKspace.jl`
- `WannierEigensystem.jl`
- `WannierWsvecIO.jl`
- `WannierWsvecGenerate.jl`
- `OrbitalProjection.jl`
- `AtomicOam.jl`
- `Win90Basis.jl`
- `HrFormat.jl`

외부 의존성:

- Julia `LinearAlgebra`
- Julia `Printf`
- Julia `TOML`
- `Spglib`

## 최소 코드 예

Lattice-reading 예:

```julia
include("spectrum/shared/CrystalCell.jl")
include("spectrum/shared/PoscarIO.jl")
include("spectrum/shared/LatticeIO.jl")

lattice = Main.LatticeIO.read_lattice("POSCAR")
println(lattice.source)
```

Projection 예:

```julia
include("spectrum/band/CLI.jl")

group = Main.BandCLI.BandCore.OrbitalProjection.ProjectionGroup("A", [1, 2], "blue")
spec = Main.BandCLI.BandCore.OrbitalProjection.ProjectionSpec([group], 4)
```

`InputParsing.jl`은 application entrypoint가 아니라 모듈별 `InputIO.jl`
파일에서 재사용하도록 만들어져 있습니다.

## 실행과 테스트

`shared`는 spectrum module check를 통해 간접적으로 테스트됩니다.

유용한 직접 확인:

```bash
julia --project=spectrum spectrum/main.jl band --help
julia --project=spectrum --startup-file=no spectrum/test/runtests.jl
```

## 흔한 실수

- `spectrum`은 구조를 내부에서 re-base하지 않습니다. HR `R` vector, reduced `KPOINTS`, manual center는 같은 input-cell basis로 맞춰야 합니다.
- fractional coordinate는 기준 cell을 함께 알 때만 의미가 있습니다.
- `band` run에서 `run.wsvec`는 명시적 opt-in입니다. 입력 옆에
  `wannier90_wsvec.dat`가 있어도 자동으로 켜지지 않습니다.
- `wsvec` validation은 numerical work 시작 전에 구조적 사용 가능성을
  확인하므로 incomplete table은 early failure가 납니다.
- `OrbitalProjection`은 Julia와 Wannier90-style orbital order에 맞춘 1-based
  Wannier index를 사용합니다.
- `AtomicOam`은 retained-basis projected local atomic OAM을 계산합니다.
  Orbital magnetization 구현은 아닙니다.
- `AtomicSpin`은 complete spinful `up/dn` basis pair를 요구합니다. SAM에서는
  unpolarized basis metadata를 error로 막습니다.

## 읽는 순서

1. `CrystalCell.jl`
2. `PoscarIO.jl`
3. `LatticeIO.jl`
4. `InputParsing.jl`
5. `CellConventions.jl`
6. `WannierTypes.jl`
7. `SpinLayout.jl`
8. `SpinExpand.jl`
9. `WannierHrIO.jl`
10. `HrHermiticity.jl`
11. `WannierKspace.jl`
12. `WannierEigensystem.jl`
13. `WannierWsvecIO.jl`
14. `WannierWsvecGenerate.jl`
15. `OrbitalProjection.jl`
16. `AtomicAngularMomentum.jl`
17. `AtomicOam.jl`
18. `Win90Basis.jl`
19. `HrFormat.jl`
