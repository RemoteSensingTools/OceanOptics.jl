// Interactive line charts for the OceanOptics.jl documentation.
//
// Pages emit a <canvas> plus a call `OOCharts.line(id, spec)` (see
// docs/DocCharts.jl). Chart.js is loaded through the RequireJS loader that
// Documenter already ships — a plain <script> tag for the UMD build would
// collide with RequireJS's anonymous-define check.
//
// Colors follow the validated reference palette: categorical slots in fixed
// order, an ordinal blue ramp for ordered series (e.g. chlorophyll levels),
// light and dark steps selected separately. The active Documenter theme is
// detected from the page's computed background, so all bundled themes work,
// and charts restyle when the theme is toggled.
(function () {
  "use strict";

  var CHARTJS = "https://cdn.jsdelivr.net/npm/chart.js@4.4.1/dist/chart.umd.min.js";

  var PALETTE = {
    light: {
      series: ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
               "#e87ba4", "#008300", "#6250d6", "#e34948"],
      ordinal: ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281", "#0d366b"],
      reference: "#52514e"
    },
    dark: {
      series: ["#3987e5", "#d95926", "#199e70", "#c98500",
               "#d55181", "#008300", "#9085e9", "#e66767"],
      ordinal: ["#184f95", "#1c5cab", "#2a78d6", "#5598e7", "#86b6ef", "#b7d3f6"],
      reference: "#c3c2b7"
    }
  };

  var charts = [];

  function rgb(css) {
    var m = css.match(/\d+(\.\d+)?/g);
    return m ? m.slice(0, 3).map(Number) : [255, 255, 255];
  }

  function isDark() {
    var bg = rgb(getComputedStyle(document.body).backgroundColor);
    return (0.2126 * bg[0] + 0.7152 * bg[1] + 0.0722 * bg[2]) / 255 < 0.5;
  }

  function ink() {
    var c = rgb(getComputedStyle(document.body).color);
    return {
      text: "rgb(" + c.join(",") + ")",
      muted: "rgba(" + c.join(",") + ",0.65)",
      grid: "rgba(" + c.join(",") + ",0.10)",
      axis: "rgba(" + c.join(",") + ",0.30)"
    };
  }

  // Evenly spaced picks from the ordinal ramp, low → high magnitude.
  function ordinal(ramp, k, n) {
    if (n <= 1) return ramp[Math.floor(ramp.length / 2)];
    return ramp[Math.round(k * (ramp.length - 1) / (n - 1))];
  }

  function seriesColor(pal, ds, k, spec) {
    if (ds.reference) return pal.reference;
    if (spec.ordinal) return ordinal(pal.ordinal, k, spec.datasets.length);
    // `slot` arrives as null when not pinned; fall back to the series order.
    var slot = typeof ds.slot === "number" ? ds.slot : k;
    return pal.series[slot % pal.series.length];
  }

  function style(chart, spec) {
    var pal = isDark() ? PALETTE.dark : PALETTE.light;
    var c = ink();
    chart.data.datasets.forEach(function (d, k) {
      var col = seriesColor(pal, spec.datasets[k], k, spec);
      d.borderColor = col;
      d.backgroundColor = col;
      d.pointBorderColor = getComputedStyle(document.body).backgroundColor;
    });
    ["x", "y"].forEach(function (ax) {
      var s = chart.options.scales[ax];
      s.ticks.color = c.muted;
      s.title.color = c.text;
      s.grid.color = c.grid;
      s.border = { color: c.axis };
    });
    chart.options.plugins.legend.labels.color = c.text;
    chart.options.plugins.title.color = c.text;
  }

  function fmt(v) {
    if (v === null || v === undefined || !isFinite(v)) return "–";
    var a = Math.abs(v);
    return (a !== 0 && (a < 1e-3 || a >= 1e5)) ? v.toExponential(4) : v.toPrecision(5);
  }

  function build(Chart, id, spec) {
    var canvas = document.getElementById(id);
    if (!canvas) return;
    var sameX = spec.datasets.every(function (d) {
      return d.x.length === spec.datasets[0].x.length;
    });
    var datasets = spec.datasets.map(function (d) {
      return {
        label: d.label,
        data: d.x.map(function (x, i) { return { x: x, y: d.y[i] }; }),
        showLine: !d.points,
        borderWidth: 2,
        borderDash: d.dash ? [6, 4] : (d.dot ? [2, 3] : []),
        pointRadius: d.points ? 4 : 0,
        pointHoverRadius: d.points ? 6 : 4,
        pointBorderWidth: d.points ? 2 : 0,
        spanGaps: false,
        tension: 0
      };
    });
    var scale = function (title, log, min, max) {
      return {
        type: log ? "logarithmic" : "linear",
        title: { display: !!title, text: title },
        min: min === null ? undefined : min,
        max: max === null ? undefined : max,
        ticks: { maxTicksLimit: 8 },
        grid: {}
      };
    };
    // Default the x axis to the data extent rather than padded "nice" ticks.
    var xs = [].concat.apply([], spec.datasets.map(function (d) { return d.x; }));
    var xmin = spec.xmin !== null ? spec.xmin : Math.min.apply(null, xs);
    var xmax = spec.xmax !== null ? spec.xmax : Math.max.apply(null, xs);
    var chart = new Chart(canvas, {
      type: "scatter",
      data: { datasets: datasets },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        animation: false,
        interaction: { mode: sameX ? "index" : "nearest", axis: "x", intersect: false },
        scales: {
          x: scale(spec.xlabel, spec.xlog, xmin, xmax),
          y: scale(spec.ylabel, spec.ylog, spec.ymin, spec.ymax)
        },
        plugins: {
          title: { display: !!spec.title, text: spec.title },
          legend: {
            display: spec.datasets.length > 1,
            position: "bottom",
            labels: { usePointStyle: true, boxHeight: 6 }
          },
          tooltip: {
            callbacks: {
              title: function (items) {
                return items.length ? (spec.xname || "x") + " = " + fmt(items[0].parsed.x) : "";
              },
              label: function (item) {
                return " " + item.dataset.label + ": " + fmt(item.parsed.y);
              }
            }
          }
        }
      }
    });
    charts.push({ chart: chart, spec: spec });
    style(chart, spec);
    chart.update();
  }

  function restyle() {
    charts.forEach(function (c) { style(c.chart, c.spec); c.chart.update(); });
  }

  new MutationObserver(restyle).observe(document.documentElement,
                                        { attributes: true, attributeFilter: ["class"] });
  if (window.matchMedia) {
    var mq = window.matchMedia("(prefers-color-scheme: dark)");
    if (mq.addEventListener) mq.addEventListener("change", restyle);
  }

  window.OOCharts = {
    line: function (id, spec) {
      require([CHARTJS], function (Chart) { build(Chart, id, spec); });
    }
  };
})();
