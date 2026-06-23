module PlotFlux3D

using Printf

using ..Model: FluxBasisEntry, FluxEdge, PlotConfig

export write_flux_html

# Don't fix this size
const FLUX_ARROW_SIZEREF = 1

function _js_string(text::AbstractString)
    escaped = replace(String(text), "\\" => "\\\\", "\"" => "\\\"", "\n" => "\\n")
    return "\"$escaped\""
end

_js_number(x::Real) = isfinite(Float64(x)) ? string(Float64(x)) : "null"
_js_vector(xs) = "[" * join((_js_number(x) for x in xs), ",") * "]"
_js_string_vector(xs) = "[" * join((_js_string(x) for x in xs), ",") * "]"

function _trace_json(; kwargs...)
    parts = String[]
    for (key, value) in kwargs
        push!(parts, string("\"", key, "\":", value))
    end
    return "{" * join(parts, ",") * "}"
end

function _unique_sites(basis::Vector{FluxBasisEntry})
    sites = NamedTuple[]
    seen = Set{Tuple{String, NTuple{3, Float64}}}()
    for entry in basis
        key = (entry.site_label, entry.center_frac)
        key in seen && continue
        push!(seen, key)
        push!(sites, (label=entry.site_label, species=entry.species, center=entry.center_frac))
    end
    return sites
end

function _translations(bounds::NTuple{3, Int})
    return [
        (rx, ry, rz)
        for rx in -bounds[1]:bounds[1]
        for ry in -bounds[2]:bounds[2]
        for rz in -bounds[3]:bounds[3]
    ]
end

function _translated_frac(center::NTuple{3, Float64}, R::NTuple{3, Int})
    return (
        center[1] + Float64(R[1]),
        center[2] + Float64(R[2]),
        center[3] + Float64(R[3]),
    )
end

_translation_label(R::NTuple{3, Int}) = "[$(R[1]),$(R[2]),$(R[3])]"

function _site_trace(lattice, basis::Vector{FluxBasisEntry}, supercell_bounds::NTuple{3, Int})
    xs, ys, zs = Float64[], Float64[], Float64[]
    labels, hover = String[], String[]
    for site in _unique_sites(basis)
        for R in _translations(supercell_bounds)
            frac = _translated_frac(site.center, R)
            pos = lattice * collect(frac)
            central = R == (0, 0, 0)
            image_label = string(site.label, " ", _translation_label(R))
            push!(xs, pos[1]); push!(ys, pos[2]); push!(zs, pos[3])
            push!(labels, central ? site.label : "")
            push!(hover, string(image_label, "<br>species=", site.species))
        end
    end
    return _trace_json(
        type=_js_string("scatter3d"),
        mode=_js_string("markers+text"),
        name=_js_string("supercell sites"),
        x=_js_vector(xs),
        y=_js_vector(ys),
        z=_js_vector(zs),
        text=_js_string_vector(labels),
        hovertext=_js_string_vector(hover),
        textposition=_js_string("top center"),
        # size and color are fixed to avoid distracting from the flux edges
        marker="{\"size\":9,\"color\":\"#1f77b4\",\"opacity\":0.75}",
        hoverinfo=_js_string("text"),
    )
end

function _cell_traces(lattice)
    traces = String[]
    origin = zeros(3)
    colors = ("#444444", "#666666", "#888888")
    for i in 1:3
        vec = lattice[:, i]
        push!(traces, _trace_json(
            type=_js_string("scatter3d"),
            mode=_js_string("lines"),
            name=_js_string("a$i"),
            x=_js_vector([origin[1], vec[1]]),
            y=_js_vector([origin[2], vec[2]]),
            z=_js_vector([origin[3], vec[3]]),
            line="{\"color\":\"$(colors[i])\",\"width\":4}",
            hoverinfo=_js_string("name"),
        ))
    end
    return traces
end

function _edge_hover(edge::FluxEdge)
    value = @sprintf("% .6f%+ .6fi", real(edge.value), imag(edge.value))
    return "nn=$(edge.nn)<br>$(edge.from_label) [$(edge.from_index)] -> " *
        "$(edge.to_label) [$(edge.to_index)]<br>R=$(edge.R)<br>value=$value"
end

function _shift_frac(point::NTuple{3, Float64}, shift::NTuple{3, Int})
    return (
        point[1] + Float64(shift[1]),
        point[2] + Float64(shift[2]),
        point[3] + Float64(shift[3]),
    )
end

function _shift_cart(point::NTuple{3, Float64}, lattice, shift::NTuple{3, Int})
    delta = lattice * [Float64(shift[1]), Float64(shift[2]), Float64(shift[3])]
    return (
        point[1] + Float64(delta[1]),
        point[2] + Float64(delta[2]),
        point[3] + Float64(delta[3]),
    )
end

function _translated_edge(edge::FluxEdge, lattice, shift::NTuple{3, Int})
    return FluxEdge(
        edge.nn,
        edge.from_index,
        edge.to_index,
        edge.R,
        edge.value,
        edge.from_label,
        edge.to_label,
        _shift_frac(edge.start_frac, shift),
        _shift_frac(edge.finish_frac, shift),
        _shift_cart(edge.start_cart, lattice, shift),
        _shift_cart(edge.finish_cart, lattice, shift),
    )
end

function _in_bounds(cell::NTuple{3, Int}, bounds::NTuple{3, Int})
    return all(i -> -bounds[i] <= cell[i] <= bounds[i], 1:3)
end

function _visual_shifts(edge::FluxEdge, bounds::NTuple{3, Int})
    limit = (
        bounds[1] + abs(edge.R[1]),
        bounds[2] + abs(edge.R[2]),
        bounds[3] + abs(edge.R[3]),
    )
    shifts = NTuple{3, Int}[]
    for sx in -limit[1]:limit[1], sy in -limit[2]:limit[2], sz in -limit[3]:limit[3]
        start_cell = (sx, sy, sz)
        finish_cell = (sx + edge.R[1], sy + edge.R[2], sz + edge.R[3])
        (_in_bounds(start_cell, bounds) || _in_bounds(finish_cell, bounds)) ||
            continue
        push!(shifts, start_cell)
    end
    return shifts
end

function _visual_edges(edges::Vector{FluxEdge}, lattice, bounds::NTuple{3, Int})
    visual = FluxEdge[]
    for edge in edges
        for shift in _visual_shifts(edge, bounds)
            push!(visual, _translated_edge(edge, lattice, shift))
        end
    end
    return visual
end

_default_edge_color(edge::FluxEdge) = imag(edge.value) >= 0 ? "#d62728" : "#2ca02c"

function _edge_style(edge::FluxEdge, plot_config::Union{Nothing, PlotConfig})
    size = Float64(FLUX_ARROW_SIZEREF)
    color = _default_edge_color(edge)
    isnothing(plot_config) && return (size=size, color=color)
    for style in plot_config.arrow_styles
        style.selector == edge.from_label || continue
        return (size=style.size, color=style.color)
    end
    return (size=size, color=color)
end

function _edge_line_trace(edge::FluxEdge, index::Integer, style)
    width = max(2.0, 6.0 * abs(edge.value))
    return _trace_json(
        type=_js_string("scatter3d"),
        mode=_js_string("lines"),
        name=_js_string("flux $index"),
        x=_js_vector([edge.start_cart[1], edge.finish_cart[1]]),
        y=_js_vector([edge.start_cart[2], edge.finish_cart[2]]),
        z=_js_vector([edge.start_cart[3], edge.finish_cart[3]]),
        line="{\"color\":$(_js_string(style.color)),\"width\":$width}",
        text="[" * _js_string(_edge_hover(edge)) * "," * _js_string(_edge_hover(edge)) * "]",
        hoverinfo=_js_string("text"),
        showlegend="false",
    )
end

function _phase_direction(edge::FluxEdge)
    im = imag(edge.value)
    im < 0 && return 1.0
    im > 0 && return -1.0
    return real(edge.value) < 0 ? 1.0 : -1.0
end

function _midpoint(a::NTuple{3, Float64}, b::NTuple{3, Float64})
    return ((a[1] + b[1]) / 2, (a[2] + b[2]) / 2, (a[3] + b[3]) / 2)
end

function _edge_cone_trace(edge::FluxEdge, index::Integer, style)
    dx = edge.finish_cart[1] - edge.start_cart[1]
    dy = edge.finish_cart[2] - edge.start_cart[2]
    dz = edge.finish_cart[3] - edge.start_cart[3]
    direction = _phase_direction(edge)
    tip = _midpoint(edge.start_cart, edge.finish_cart)
    return _trace_json(
        type=_js_string("cone"),
        name=_js_string("flux direction $index"),
        x=_js_vector([tip[1]]),
        y=_js_vector([tip[2]]),
        z=_js_vector([tip[3]]),
        u=_js_vector([direction * dx]),
        v=_js_vector([direction * dy]),
        w=_js_vector([direction * dz]),
        text="[" * _js_string(_edge_hover(edge)) * "]",
        hoverinfo=_js_string("text"),
        sizemode=_js_string("absolute"),
        sizeref=_js_number(style.size),
        anchor=_js_string("tip"),
        colorscale="[[0,$(_js_string(style.color))],[1,$(_js_string(style.color))]]",
        showscale="false",
        showlegend="false",
    )
end

function _plotly_page(traces::Vector{String})
    layout = """
{"title":"Flux bonds","margin":{"l":0,"r":0,"b":0,"t":42},
"scene":{"aspectmode":"data","xaxis":{"title":"x"},"yaxis":{"title":"y"},"zaxis":{"title":"z"}}}
"""
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

function write_flux_html(
    path::AbstractString,
    lattice::AbstractMatrix{<:Real},
    basis::Vector{FluxBasisEntry},
    edges::Vector{FluxEdge},
    supercell_bounds::NTuple{3, Int}=(0, 0, 0),
    plot_config::Union{Nothing, PlotConfig}=nothing,
)
    visual_edges = _visual_edges(edges, lattice, supercell_bounds)
    traces = String[_site_trace(lattice, basis, supercell_bounds)]
    append!(traces, _cell_traces(lattice))
    for (i, edge) in enumerate(visual_edges)
        style = _edge_style(edge, plot_config)
        push!(traces, _edge_line_trace(edge, i, style))
        push!(traces, _edge_cone_trace(edge, i, style))
    end
    mkpath(dirname(abspath(path)))
    write(path, _plotly_page(traces))
    return path
end

end
