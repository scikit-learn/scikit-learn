"""
Tools for matplotlib plot exporting
"""


def ipynb_vega_init():
    """Initialize the IPython notebook display elements

    This function borrows heavily from the excellent vincent package:
    http://github.com/wrobstory/vincent
    """
    try:
        from IPython.core.display import display, HTML
    except ImportError:
        print("IPython Notebook could not be loaded.")

    require_js = """
    if (window['d3'] === undefined) {{
        require.config({{ paths: {{d3: "http://d3js.org/d3.v3.min"}} }});
        require(["d3"], function(d3) {{
          window.d3 = d3;
          {0}
        }});
    }};
    if (window['topojson'] === undefined) {{
        require.config(
            {{ paths: {{topojson: "http://d3js.org/topojson.v1.min"}} }}
            );
        require(["topojson"], function(topojson) {{
          window.topojson = topojson;
        }});
    }};
    """
    d3_geo_projection_js_url = "http://d3js.org/d3.geo.projection.v0.min.js"
    d3_layout_cloud_js_url = "http://wrobstory.github.io/d3-cloud/d3.layout.cloud.js"
    topojson_js_url = "http://d3js.org/topojson.v1.min.js"
    vega_js_url = "http://trifacta.github.com/vega/vega.js"

    dep_libs = """$.getScript("%s", function() {
        $.getScript("%s", function() {
            $.getScript("%s", function() {
                $.getScript("%s", function() {
                        $([IPython.events]).trigger("vega_loaded.vincent");
                })
            })
        })
    });""" % (
        d3_geo_projection_js_url,
        d3_layout_cloud_js_url,
        topojson_js_url,
        vega_js_url,
    )
    load_js = require_js.format(dep_libs)
    html = "<script>" + load_js + "</script>"
    display(HTML(html))
