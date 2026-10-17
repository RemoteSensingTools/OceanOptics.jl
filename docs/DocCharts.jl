"""
    DocCharts

Interactive charts for the documentation. `linechart` returns an `HTML`
object that an `@example` block can display: a Chart.js canvas (rendered by
`docs/src/assets/oo_charts.js`) followed by a collapsible table of the
plotted values, so every figure is also readable as data.
"""
module DocCharts

export Series, linechart

"""
    Series(label, x, y; points = false, dash = false, dot = false,
           reference = false, slot = nothing)

One plotted series. `points = true` draws markers without a connecting line
(reference data); `reference = true` uses the neutral reference ink instead
of a categorical color; `slot` pins a categorical color slot (0-based) so an
entity keeps its color across charts.
"""
struct Series
    label::String
    x::Vector{Float64}
    y::Vector{Float64}
    points::Bool
    dash::Bool
    dot::Bool
    reference::Bool
    slot::Union{Nothing, Int}
end

Series(label, x, y; points = false, dash = false, dot = false, reference = false,
       slot = nothing) =
    Series(string(label), collect(Float64, x), collect(Float64, y), points, dash, dot,
           reference, slot)

const _counter = Ref(0)

"""
    linechart(series; xlabel, ylabel, title = "", xname = "x", xlog = false,
              ylog = false, xmin, xmax, ymin, ymax, ordinal = false,
              height = 380, table = true, digits = 5) -> HTML

`ordinal = true` colors the series with an ordered one-hue ramp (for series
that are levels of one quantity, such as chlorophyll concentrations);
otherwise categorical slots are assigned in series order.
"""
function linechart(series::AbstractVector{Series}; xlabel = "", ylabel = "", title = "",
                   xname = "x", xlog = false, ylog = false, xmin = nothing,
                   xmax = nothing, ymin = nothing, ymax = nothing, ordinal = false,
                   height = 380, table = true, digits = 5)
    id   = "oo-chart-$(_counter[] += 1)"
    spec = Dict{String, Any}(
        "xlabel" => xlabel, "ylabel" => ylabel, "title" => title, "xname" => xname,
        "xlog" => xlog, "ylog" => ylog, "xmin" => xmin, "xmax" => xmax,
        "ymin" => ymin, "ymax" => ymax, "ordinal" => ordinal,
        "datasets" => [Dict{String, Any}("label" => s.label, "x" => s.x, "y" => s.y,
                                         "points" => s.points, "dash" => s.dash,
                                         "dot" => s.dot, "reference" => s.reference,
                                         "slot" => s.slot) for s in series])
    io = IOBuffer()
    print(io, """<div class="oo-chart" style="position:relative;height:$(height)px;margin:1em 0">""",
          """<canvas id="$id" role="img" aria-label="$(_html(title == "" ? ylabel : title))"></canvas></div>""",
          """<script>OOCharts.line("$id", """, _json(spec), ");</script>")
    table && _table(io, series, xname, digits)
    return HTML(String(take!(io)))
end

# Collapsible data table: one row per x value of the first series, one
# column per series (blank where a series has no value at that x).
function _table(io, series, xname, digits)
    xs = sort(unique(reduce(vcat, (s.x for s in series))))
    print(io, "<details><summary>Data table</summary><table><thead><tr><th>",
          _html(xname), "</th>")
    foreach(s -> print(io, "<th>", _html(s.label), "</th>"), series)
    print(io, "</tr></thead><tbody>")
    for x in xs
        print(io, "<tr><td>", _num(x, digits), "</td>")
        for s in series
            i = findfirst(==(x), s.x)
            print(io, "<td>", i === nothing ? "" : _num(s.y[i], digits), "</td>")
        end
        print(io, "</tr>")
    end
    print(io, "</tbody></table></details>")
end

_num(x, digits) = isfinite(x) ? string(round(x; sigdigits = digits)) : ""
_html(s) = replace(string(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "\"" => "&quot;")

# Minimal JSON writer for the chart spec (numbers, strings, Bools, nothing,
# vectors, string-keyed dicts). Plotted numbers are rounded to six
# significant digits to keep pages small; the data table shows `digits`.
# "</" is escaped so the payload cannot close the surrounding <script>.
_json(x::Nothing) = "null"
_json(x::Bool) = x ? "true" : "false"
_json(x::Integer) = string(x)
_json(x::AbstractFloat) = isfinite(x) ? repr(round(Float64(x); sigdigits = 6)) : "null"
_json(x::AbstractString) = "\"" * replace(x, "\\" => "\\\\", "\"" => "\\\"", "\n" => "\\n",
                                          "</" => "<\\/") * "\""
_json(x::AbstractVector) = "[" * join((_json(v) for v in x), ",") * "]"
_json(x::AbstractDict) = "{" * join((_json(string(k)) * ":" * _json(v) for (k, v) in x), ",") * "}"

end # module
