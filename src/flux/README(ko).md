# flux

언어: [English](README.md) | 한국어

`flux`는 기존 Wannier `hr.dat`에 directed complex flux term을 추가하고, seed
bond만 Plotly 3D HTML로 보여주는 `spectrum/` 내부 모듈입니다. 구현은
`spectrum/shared` 안의 지원 모듈을 사용합니다.

## 사용법

```bash
julia spectrum/main.jl flux --input input.flux.toml --output flux_hr.dat --html flux.html --diagnostic flux_diagnostic.tsv
```

입력 예:

```toml
[flux.run]
hr = "wannier90_hr.dat"
win = "wannier90.win"
# poscar = "POSCAR"  # win이 없을 때 사용

[flux.geometry]
search_bounds = [1, 1, 0]
distance_tol = 1.0e-8

[flux.plot]
interactive = true
cell_bounds = [1, 1, 0]
arrow_styles = [
  ["C1:s", 1.0, "#2ca02c"],
  ["C2:s", 1.0, "#d62728"],
]

# POSCAR fallback에서 species별 orbital 수가 다를 때 선택적으로 지정합니다.
[flux.basis]
orbitals_per_atom = [["Mo", 5], ["S", 3]]
# orbitals_per_atom = [5, 3]  # POSCAR species-group 순서; 중복 label에 유용

[flux.diagnostic]
continuity = true
continuity_tol = 1.0e-10
plaquettes = [
  ["C1-triangle", [
    ["C1", [0, 0, 0]],
    ["C1", [1, 1, 0]],
    ["C1", [1, 0, 0]],
  ]],
]

[[flux.terms]]
term = [
  [2, ["C1", "C1"], [1, 1, 0], [0.0, -0.02]],
  [2, ["C1", "C1"], [1, 0, 0], [0.0, 0.02]],
  [2, ["C1", "C1"], [0, 1, 0], [0.0, 0.02]],
]
```

`nn`은 distance-shell index이며 물리적 neighbor label과 항상 같지는 않습니다.
숫자 endpoint는 1-based Wannier index입니다. 문자열 endpoint는 먼저 `Mo1` 같은
`site_label`로 해석하고, 없으면 `Mo` 같은 species 전체로 확장합니다.
각 term row는 `[nn, [from, to], [rx, ry, rz], [re, im]]` 형식이므로 seed cell
translation `R`을 직접 지정합니다. broad selector가 Hermitian-equivalent seed를
같이 잡으면 하나만 남기고, 서로 다른 row가 같은 physical pair를 다시 지정하면
error를 냅니다.

출력 경로는 TOML key가 아니라 CLI option으로 지정합니다. 출력 HR은
`ndegen = 1`인 normalized 형식이며, HTML에는 Hermitian partner가 아닌 seed flux
bond와 그 incoming periodic image를 표시합니다. row는 Wannier matrix element
`<from,0|H|to,R>`에 대응하므로 양의 phase 방향은 `to@R -> from@0`입니다. 방향
cone은 bond midpoint에 두고, 3D view에서 bond 방향이 잘 보이도록 큰 arrowhead로
그립니다. imaginary phase가 음수이면 visual arrow를 반대로, 즉 `from@0 -> to@R`로
그립니다.
`[flux.plot] cell_bounds`는 표시할 supercell image 범위를 정하며, atom marker와
반복된 periodic seed arrow 모두에 적용됩니다.
생략하면 `[flux.geometry] search_bounds`를 사용합니다.
`arrow_styles` row는 `[from_atom_orbital_label, arrow_size, color]` 형식으로
directed seed arrow와 그 incoming periodic image에 적용됩니다.

`[flux.diagnostic]`는 기본적으로 `outputs/diagnostic/` 아래 TSV report를 씁니다.
경로를 바꾸려면 `--diagnostic path.tsv`, 쓰지 않으려면 `--no-diagnostic`을 사용합니다.
plaquette diagnostic은 matched seed flux edge의 `imag(value)`만 합산하며, 최종
Hamiltonian Wilson-loop phase를 계산하지 않습니다. `continuity = true`는 arrow와 같은
imaginary-sign 방향 convention으로 orbital별 flow in/out residual도 보고합니다.

`win`이 없으면 `poscar`를 지정할 수 있습니다. POSCAR fallback은 projection
metadata가 없기 때문에 HR index를 POSCAR atom 순서에 따라 배정합니다.
`[flux.basis]`가 없으면 `num_wann`이 atom 수로 균등 분배될 때만 허용합니다.
중복 species row는 count가 같으면 허용합니다. POSCAR species group별 count가
다르면 POSCAR group 순서의 integer vector 형식을 쓰면 됩니다. 원소명이 없는
VASP4-style POSCAR는 `Type1`, `Type2` synthetic species를 쓰고 site label은
`Type1_1`처럼 생성합니다.
