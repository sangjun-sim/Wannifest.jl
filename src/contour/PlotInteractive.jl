module PlotInteractive

using ..Model: EnergySurfaceResult, RunConfig
using ..Output

export open_interactive

const TRACE_COLORS = (
    "#1f77b4",
    "#d62728",
    "#2ca02c",
    "#9467bd",
    "#ff7f0e",
    "#17becf",
)

_axis_label(axis::Int) = axis == 1 ? "kx" : axis == 2 ? "ky" : "kz"
_band_label(result::EnergySurfaceResult, slot::Int) = "band $(result.bands[slot])"
_js_number(x::Real) = isfinite(Float64(x)) ? string(Float64(x)) : "null"
_js_vector(xs) = "[" * join((_js_number(x) for x in xs), ",") * "]"

function _xy_labels(config::RunConfig)
    return _axis_label(config.plane.x_axis), _axis_label(config.plane.y_axis)
end

function _js_string(text::AbstractString)
    escaped = replace(String(text), "\\" => "\\\\", "\"" => "\\\"", "\n" => "\\n")
    return "\"$escaped\""
end

function _js_matrix(z)
    rows = String[]
    for iy in axes(z, 1)
        push!(rows, _js_vector(@view z[iy, :]))
    end
    return "[" * join(rows, ",") * "]"
end

function _plotly_colorscale(config::RunConfig)
    lower = lowercase(config.plot.colormap)
    lower == "viridis" && return "Viridis"
    lower == "plasma" && return "Plasma"
    lower == "inferno" && return "Inferno"
    lower == "magma" && return "Magma"
    return config.plot.colormap
end

function _plotly_page(traces::Vector{String}, layout::String)
    return """
<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
  <style>html, body, #plot { width: 100%; height: 100%; margin: 0; }</style>
</head>
<body>
  <div id="plot"></div>
  <script>
    const data = [$(join(traces, ","))];
    const layout = $layout;
    Plotly.newPlot("plot", data, layout, {
      responsive: true,
      scrollZoom: true,
      displaylogo: false
    });
  </script>
</body>
</html>
"""
end

function _write_interactive_html(path::AbstractString, traces::Vector{String}, layout::String)
    mkpath(dirname(path))
    write(path, _plotly_page(traces, layout))
    return path
end

function _surface_html(result::EnergySurfaceResult, config::RunConfig)
    traces = String[]
    for band_slot in eachindex(result.bands)
        push!(traces, """
{"type":"surface","name":$(_js_string(_band_label(result, band_slot))),
"x":$(_js_vector(result.x_axis)),"y":$(_js_vector(result.y_axis)),
"z":$(_js_matrix(result.energies[:, :, band_slot])),
"opacity":$(band_slot == 1 ? 0.95 : 0.68),
"showscale":$(band_slot == length(result.bands)),
"colorscale":$(_js_string(_plotly_colorscale(config)))}
""")
    end
    xlabel, ylabel = _xy_labels(config)
    layout = """
{"title":"Energy surface","margin":{"l":0,"r":0,"b":0,"t":42},
"scene":{"xaxis":{"title":$(_js_string(xlabel))},
"yaxis":{"title":$(_js_string(ylabel))},
"zaxis":{"title":"E (eV)","range":$(_js_vector(config.plot.energy_range))}}}
"""
    return _write_interactive_html(Output.default_interactive_path(config, :surface), traces, layout)
end

function _contour_html(result::EnergySurfaceResult, config::RunConfig)
    traces = String[]
    for band_slot in eachindex(result.bands)
        color = TRACE_COLORS[mod1(band_slot, length(TRACE_COLORS))]
        push!(traces, """
{"type":"contour","name":$(_js_string(_band_label(result, band_slot))),
"x":$(_js_vector(result.x_axis)),"y":$(_js_vector(result.y_axis)),
"z":$(_js_matrix(result.energies[:, :, band_slot])),
"ncontours":$(config.plot.contour_levels),
"contours":{"coloring":"lines","showlabels":true},
"line":{"color":$(_js_string(color)),"width":2},"showscale":false}
""")
    end
    xlabel, ylabel = _xy_labels(config)
    layout = """
{"title":"Energy contours","xaxis":{"title":$(_js_string(xlabel))},
"yaxis":{"title":$(_js_string(ylabel))},"legend":{"orientation":"h"},
"margin":{"l":64,"r":24,"b":56,"t":48}}
"""
    return _write_interactive_html(Output.default_interactive_path(config, :contour), traces, layout)
end

function _grid_layout(n::Int)
    cols = ceil(Int, sqrt(n))
    rows = ceil(Int, n / cols)
    return rows, cols
end

function _axis_name(prefix::AbstractString, slot::Int)
    slot == 1 && return prefix
    return string(prefix, slot)
end

function _domain(index::Int, rows::Int, cols::Int)
    row = fld(index - 1, cols) + 1
    col = mod1(index, cols)
    gap = 0.04
    return (
        (col - 1) / cols + gap / 2,
        col / cols - gap / 2,
        1 - row / rows + gap / 2,
        1 - (row - 1) / rows - gap / 2,
    )
end

function _heatmap_html(result::EnergySurfaceResult, config::RunConfig)
    rows, cols = _grid_layout(length(result.bands))
    traces, axes, annotations = String[], String[], String[]
    xlabel, ylabel = _xy_labels(config)
    for band_slot in eachindex(result.bands)
        x0, x1, y0, y1 = _domain(band_slot, rows, cols)
        xaxis, yaxis = _axis_name("x", band_slot), _axis_name("y", band_slot)
        xlayout, ylayout = _axis_name("xaxis", band_slot), _axis_name("yaxis", band_slot)
        push!(traces, """
{"type":"heatmap","name":$(_js_string(_band_label(result, band_slot))),
"x":$(_js_vector(result.x_axis)),"y":$(_js_vector(result.y_axis)),
"z":$(_js_matrix(result.energies[:, :, band_slot])),
"xaxis":$(_js_string(xaxis)),"yaxis":$(_js_string(yaxis)),
"zmin":$(config.plot.energy_range[1]),"zmax":$(config.plot.energy_range[2]),
"colorscale":$(_js_string(_plotly_colorscale(config))),
"showscale":$(band_slot == length(result.bands))}
""")
        push!(axes, "\"$xlayout\":{\"domain\":[$x0,$x1],\"title\":$(_js_string(xlabel))}")
        push!(axes, "\"$ylayout\":{\"domain\":[$y0,$y1],\"title\":$(_js_string(ylabel))}")
        push!(annotations, """
{"text":$(_js_string(_band_label(result, band_slot))),"x":$((x0 + x1) / 2),
"y":$y1,"xref":"paper","yref":"paper","showarrow":false,"yshift":12}
""")
    end
    layout = """
{"title":"Energy heatmaps","margin":{"l":64,"r":42,"b":56,"t":64},
$(join(axes, ",\n")),"annotations":[$(join(annotations, ","))]}
"""
    return _write_interactive_html(Output.default_interactive_path(config, :heatmap), traces, layout)
end

function _open_in_browser(path::AbstractString)
    absolute = abspath(path)
    cmd = Sys.isapple() ? `open $absolute` :
        Sys.iswindows() ? Cmd(["cmd", "/c", "start", "", absolute]) :
        `xdg-open $absolute`
    try
        run(cmd)
    catch err
        @warn "Could not open interactive HTML in a browser" path=absolute exception=(err, catch_backtrace())
    end
    return nothing
end

function open_interactive(
    result::EnergySurfaceResult,
    config::RunConfig;
    open_browser::Bool=true,
    opener=_open_in_browser,
)
    paths = String[]
    config.plot.mode in (:surface, :both) && push!(paths, _surface_html(result, config))
    config.plot.mode in (:contour, :both) && push!(paths, _contour_html(result, config))
    config.plot.mode == :heatmap && push!(paths, _heatmap_html(result, config))
    open_browser && foreach(opener, paths)
    return paths
end

end
