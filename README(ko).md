# Wannifest.jl

Wannifest는 선택된 spectrum workflow를 복사해 둔 Julia package입니다.

- `band`: band structure, DOS, projection, OAM, SAM workflow.
- `contour`: 2D k-plane energy surface data와 plot.
- `flux`: Wannier Hamiltonian에 directed flux term 추가.
- `supercell`: POSCAR/CONTCAR structural supercell builder.
- `superham`: supercell Hamiltonian builder.

원본 저장소는 그대로 둡니다. 이 package는 복사된 command가 필요로 하는
shared helper layer도 의도적으로 함께 복사합니다.

## 사용법

```bash
julia --project=/Users/sangjun/Wannifest.jl /Users/sangjun/Wannifest.jl/main.jl --help
julia --project=/Users/sangjun/Wannifest.jl /Users/sangjun/Wannifest.jl/main.jl band --input /Users/sangjun/Wannifest.jl/examples/flake/graphene/input.band.toml --no-plot
julia --project=/Users/sangjun/Wannifest.jl /Users/sangjun/Wannifest.jl/main.jl contour --input /Users/sangjun/Wannifest.jl/examples/contour/graphene/input.toml --no-plot
julia --project=/Users/sangjun/Wannifest.jl /Users/sangjun/Wannifest.jl/main.jl superham --input /Users/sangjun/Wannifest.jl/examples/flake/graphene/input.superham.toml
julia --project=/Users/sangjun/Wannifest.jl /Users/sangjun/Wannifest.jl/main.jl supercell 2 0 0, 0 2 0, 0 0 1
julia --project=/Users/sangjun/Wannifest.jl /Users/sangjun/Wannifest.jl/main.jl flux --input /Users/sangjun/Wannifest.jl/examples/flux/graphene/input.toml --no-html --no-diagnostic
```

## 포함하지 않는 것

이 복사본에는 flake, SOC, utility command가 없습니다.

## 테스트

```bash
julia --project=/Users/sangjun/Wannifest.jl --startup-file=no /Users/sangjun/Wannifest.jl/test/runtests.jl
```
