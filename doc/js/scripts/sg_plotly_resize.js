// Related to https://github.com/scikit-learn/scikit-learn/issues/30279
// There an interaction between plotly and bootstrap/pydata-sphinx-theme
// that causes plotly figures to not detect the right-hand sidebar width

function resizePlotlyGraphs() {
    const plotlyDivs = document.getElementsByClassName("plotly-graph-div");

    for (const div of plotlyDivs) {
        Plotly.Plots.resize(div);
    }
}

window.addEventListener("resize", resizePlotlyGraphs);
document.addEventListener("DOMContentLoaded", resizePlotlyGraphs);
