# spectrum contour

`spectrum/contour`는 reduced-coordinate k-plane을 2차원으로 샘플링하고
`(kx, ky, kz, band, energy)` 데이터와 선택적 PNG plot을 씁니다.

저장소 루트에서 실행:

```bash
julia spectrum/main.jl contour --input path/to/input.contour.toml
```

plot 생성을 건너뛰기:

```bash
julia spectrum/main.jl contour --input path/to/input.contour.toml --no-plot
```

출력 디렉터리 변경:

```bash
julia spectrum/main.jl contour --input path/to/input.contour.toml --output-dir outputs/custom
```

## 입력

예시:

```toml
[contour.run]
hr = "wannier90_hr.dat"
structure = ""
wsvec = ""
hermiticity_tol = 1.0e-8
verbose = true

[contour.plane]
axes = ["kx", "ky"]
fixed_axis = "kz"
fixed_value = 0.0
range_x = [-0.5, 0.5]
range_y = [-0.5, 0.5]
mesh = "101x101"

[contour.energy]
shift = 0.0
bands = [1]

[contour.spin]
layout = "qe"

[contour.plot]
mode = "both"
interactive = false
size = [900, 700]
energy_range = [-3.0, 3.0]
colormap = "viridis"
contour_levels = 40
```

입력 경로는 TOML 파일 위치를 기준으로 해석합니다. 출력 파일 경로는 TOML schema에
넣지 않으며, 기본 디렉터리를 옮겨야 할 때만 `--output-dir`를 사용합니다.

`[contour.spin].layout`은 source HR/wsvec ordering(`qe` 또는 `vasp544`)만
뜻합니다. 현재 contour 모듈은 sorted band index surface를 시각화하며 spin
expectation map은 계산하지 않습니다.

## 출력

기본 출력은 입력 TOML 디렉터리 기준으로 작성됩니다.

- `outputs/data/contour_energy_surface.dat`
- `outputs/plots/contour_surface.png`
- `outputs/plots/contour_contour.png`
- `outputs/plots/contour_heatmap.png`

`interactive = true`이면 대응되는 Plotly HTML 파일도 저장하고 기본 브라우저로
엽니다.

- `outputs/plots/contour_surface.html`
- `outputs/plots/contour_contour.html`
- `outputs/plots/contour_heatmap.html`

브라우저 plot은 zoom, pan, hover를 지원하고 surface mode에서는 3D rotation도
가능합니다. `--no-plot`은 PNG 저장과 interactive HTML 열기를 모두 끕니다.
선택한 모든 band는 plot mode별 한 figure 안에 함께 그립니다.
