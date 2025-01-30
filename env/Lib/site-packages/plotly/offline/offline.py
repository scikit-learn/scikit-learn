""" Plotly Offline
    A module to use Plotly's graphing library with Python
    without connecting to a public or private plotly enterprise
    server.
"""
import os
import warnings
import pkgutil

from plotly.optional_imports import get_module
from plotly import tools
from ._plotlyjs_version import __plotlyjs_version__


__IMAGE_FORMATS = ["jpeg", "png", "webp", "svg"]


def download_plotlyjs(download_url):
    warnings.warn(
        """
        `download_plotlyjs` is deprecated and will be removed in the
        next release. plotly.js is shipped with this module, it is no
        longer necessary to download this bundle separately.
    """,
        DeprecationWarning,
    )


def get_plotlyjs_version():
    """
    Returns the version of plotly.js that is bundled with plotly.py.

    Returns
    -------
    str
        Plotly.js version string
    """
    return __plotlyjs_version__


def get_plotlyjs():
    """
    Return the contents of the minified plotly.js library as a string.

    This may be useful when building standalone HTML reports.

    Returns
    -------
    str
        Contents of the minified plotly.js library as a string

    Examples
    --------
    Here is an example of creating a standalone HTML report that contains
    two plotly figures, each in their own div.  The include_plotlyjs argument
    is set to False when creating the divs so that we don't include multiple
    copies of the plotly.js library in the output.  Instead, a single copy
    of plotly.js is included in a script tag in the html head element.

    >>> import plotly.graph_objs as go
    >>> from plotly.offline import plot, get_plotlyjs
    >>> fig1 = go.Figure(data=[{'type': 'bar', 'y': [1, 3, 2]}],
    ...                 layout={'height': 400})
    >>> fig2 = go.Figure(data=[{'type': 'scatter', 'y': [1, 3, 2]}],
    ...                  layout={'height': 400})
    >>> div1 = plot(fig1, output_type='div', include_plotlyjs=False)
    >>> div2 = plot(fig2, output_type='div', include_plotlyjs=False)

    >>> html = '''
    ... <html>
    ...     <head>
    ...         <script type="text/javascript">{plotlyjs}</script>
    ...     </head>
    ...     <body>
    ...        {div1}
    ...        {div2}
    ...     </body>
    ... </html>
    ... '''.format(plotlyjs=get_plotlyjs(), div1=div1, div2=div2)

    >>> with open('multi_plot.html', 'w') as f:
    ...      f.write(html) # doctest: +SKIP
    """
    path = os.path.join("package_data", "plotly.min.js")
    plotlyjs = pkgutil.get_data("plotly", path).decode("utf-8")
    return plotlyjs


def _build_resize_script(plotdivid, plotly_root="Plotly"):
    resize_script = (
        '<script type="text/javascript">'
        'window.addEventListener("resize", function(){{'
        'if (document.getElementById("{id}")) {{'
        '{plotly_root}.Plots.resize(document.getElementById("{id}"));'
        "}};}})"
        "</script>"
    ).format(plotly_root=plotly_root, id=plotdivid)
    return resize_script


def _build_mathjax_script(url):
    return '<script src="{url}?config=TeX-AMS-MML_SVG"></script>'.format(url=url)


def _get_jconfig(config=None):

    configkeys = (
        "staticPlot",
        "plotlyServerURL",
        "editable",
        "edits",
        "autosizable",
        "responsive",
        "queueLength",
        "fillFrame",
        "frameMargins",
        "scrollZoom",
        "doubleClick",
        "showTips",
        "showAxisDragHandles",
        "showAxisRangeEntryBoxes",
        "showLink",
        "sendData",
        "showSendToCloud",
        "linkText",
        "showSources",
        "displayModeBar",
        "modeBarButtonsToRemove",
        "modeBarButtonsToAdd",
        "modeBarButtons",
        "toImageButtonOptions",
        "displaylogo",
        "watermark",
        "plotGlPixelRatio",
        "setBackground",
        "topojsonURL",
        "mapboxAccessToken",
        "logging",
        "globalTransforms",
        "locale",
        "locales",
        "doubleClickDelay",
    )

    if config and isinstance(config, dict):
        # Warn user on unrecognized config options.  We make this a warning
        # rather than an error since we don't have code generation logic in
        # place yet to guarantee that the config options in plotly.py are up
        # to date
        bad_config = [k for k in config if k not in configkeys]
        if bad_config:
            warnings.warn(
                """
Unrecognized config options supplied: {bad_config}""".format(
                    bad_config=bad_config
                )
            )

        clean_config = config
    else:
        clean_config = {}

    plotly_platform_url = tools.get_config_plotly_server_url()

    if not clean_config.get("plotlyServerURL", None):
        clean_config["plotlyServerURL"] = plotly_platform_url

    if (
        plotly_platform_url != "https://plot.ly"
        and clean_config.get("linkText", None) == "Export to plot.ly"
    ):
        link_domain = plotly_platform_url.replace("https://", "").replace("http://", "")
        link_text = clean_config["linkText"].replace("plot.ly", link_domain)
        clean_config["linkText"] = link_text

    return clean_config


# Build script to set global PlotlyConfig object. This must execute before
# plotly.js is loaded.
_window_plotly_config = """\
<script type="text/javascript">\
window.PlotlyConfig = {MathJaxConfig: 'local'};\
</script>"""

_mathjax_config = """\
<script type="text/javascript">\
if (window.MathJax && window.MathJax.Hub && window.MathJax.Hub.Config) {window.MathJax.Hub.Config({SVG: {font: "STIX-Web"}});}\
</script>"""


def get_image_download_script(caller):
    """
    This function will return a script that will download an image of a Plotly
    plot.

    Keyword Arguments:
    caller ('plot', 'iplot') -- specifies which function made the call for the
        download script. If `iplot`, then an extra condition is added into the
        download script to ensure that download prompts aren't initiated on
        page reloads.
    """

    if caller == "iplot":
        check_start = "if(document.readyState == 'complete') {{"
        check_end = "}}"
    elif caller == "plot":
        check_start = ""
        check_end = ""
    else:
        raise ValueError("caller should only be one of `iplot` or `plot`")

    return (
        "function downloadimage(format, height, width,"
        " filename) {{"
        "var p = document.getElementById('{{plot_id}}');"
        "Plotly.downloadImage(p, {{format: format, height: height, "
        "width: width, filename: filename}});}};"
        + check_start
        + "downloadimage('{format}', {height}, {width}, "
        "'{filename}');" + check_end
    )


def build_save_image_post_script(
    image, image_filename, image_height, image_width, caller
):
    if image:
        if image not in __IMAGE_FORMATS:
            raise ValueError(
                "The image parameter must be one of the "
                "following: {}".format(__IMAGE_FORMATS)
            )

        script = get_image_download_script(caller)
        post_script = script.format(
            format=image,
            width=image_width,
            height=image_height,
            filename=image_filename,
        )
    else:
        post_script = None

    return post_script


def init_notebook_mode(connected=False):
    """
    Initialize plotly.js in the browser if it hasn't been loaded into the DOM
    yet. This is an idempotent method and can and should be called from any
    offline methods that require plotly.js to be loaded into the notebook dom.

    Keyword arguments:

    connected (default=False) -- If True, the plotly.js library will be loaded
    from an online CDN. If False, the plotly.js library will be loaded locally
    from the plotly python package

    Use `connected=True` if you want your notebooks to have smaller file sizes.
    In the case where `connected=False`, the entirety of the plotly.js library
    will be loaded into the notebook, which will result in a file-size increase
    of a couple megabytes. Additionally, because the library will be downloaded
    from the web, you and your viewers must be connected to the internet to be
    able to view charts within this notebook.

    Use `connected=False` if you want you and your collaborators to be able to
    create and view these charts regardless of the availability of an internet
    connection. This is the default option since it is the most predictable.
    Note that under this setting the library will be included inline inside
    your notebook, resulting in much larger notebook sizes compared to the case
    where `connected=True`.
    """
    import plotly.io as pio

    ipython = get_module("IPython")
    if not ipython:
        raise ImportError("`iplot` can only run inside an IPython Notebook.")

    if connected:
        pio.renderers.default = "plotly_mimetype+notebook_connected"
    else:
        pio.renderers.default = "plotly_mimetype+notebook"

    # Trigger immediate activation of notebook. This way the plotly.js
    # library reference is available to the notebook immediately
    pio.renderers._activate_pending_renderers()


def iplot(
    figure_or_data,
    show_link=False,
    link_text="Export to plot.ly",
    validate=True,
    image=None,
    filename="plot_image",
    image_width=800,
    image_height=600,
    config=None,
    auto_play=True,
    animation_opts=None,
):
    """
    Draw plotly graphs inside an IPython or Jupyter notebook

    figure_or_data -- a plotly.graph_objs.Figure or plotly.graph_objs.Data or
                      dict or list that describes a Plotly graph.
                      See https://plot.ly/python/ for examples of
                      graph descriptions.

    Keyword arguments:
    show_link (default=False) -- display a link in the bottom-right corner of
                                of the chart that will export the chart to
                                Plotly Cloud or Plotly Enterprise
    link_text (default='Export to plot.ly') -- the text of export link
    validate (default=True) -- validate that all of the keys in the figure
                               are valid? omit if your version of plotly.js
                               has become outdated with your version of
                               graph_reference.json or if you need to include
                               extra, unnecessary keys in your figure.
    image (default=None |'png' |'jpeg' |'svg' |'webp') -- This parameter sets
        the format of the image to be downloaded, if we choose to download an
        image. This parameter has a default value of None indicating that no
        image should be downloaded. Please note: for higher resolution images
        and more export options, consider using plotly.io.write_image. See
        https://plot.ly/python/static-image-export/ for more details.
    filename (default='plot') -- Sets the name of the file your image
        will be saved to. The extension should not be included.
    image_height (default=600) -- Specifies the height of the image in `px`.
    image_width (default=800) -- Specifies the width of the image in `px`.
    config (default=None) -- Plot view options dictionary. Keyword arguments
        `show_link` and `link_text` set the associated options in this
        dictionary if it doesn't contain them already.
    auto_play (default=True) -- Whether to automatically start the animation
        sequence on page load, if the figure contains frames. Has no effect if
        the figure does not contain frames.
    animation_opts (default=None) -- Dict of custom animation parameters that
        are used for the automatically started animation on page load. This
        dict is passed to the function Plotly.animate in Plotly.js. See
        https://github.com/plotly/plotly.js/blob/master/src/plots/animation_attributes.js
        for available options. Has no effect if the figure
        does not contain frames, or auto_play is False.

    Example:
    ```
    from plotly.offline import init_notebook_mode, iplot
    init_notebook_mode()
    iplot([{'x': [1, 2, 3], 'y': [5, 2, 7]}])
    # We can also download an image of the plot by setting the image to the
    format you want. e.g. `image='png'`
    iplot([{'x': [1, 2, 3], 'y': [5, 2, 7]}], image='png')
    ```

    animation_opts Example:
    ```
    from plotly.offline import iplot
    figure = {'data': [{'x': [0, 1], 'y': [0, 1]}],
              'layout': {'xaxis': {'range': [0, 5], 'autorange': False},
                         'yaxis': {'range': [0, 5], 'autorange': False},
                         'title': 'Start Title'},
              'frames': [{'data': [{'x': [1, 2], 'y': [1, 2]}]},
                         {'data': [{'x': [1, 4], 'y': [1, 4]}]},
                         {'data': [{'x': [3, 4], 'y': [3, 4]}],
                          'layout': {'title': 'End Title'}}]}
    iplot(figure, animation_opts={'frame': {'duration': 1}})
    ```
    """
    import plotly.io as pio

    ipython = get_module("IPython")
    if not ipython:
        raise ImportError("`iplot` can only run inside an IPython Notebook.")

    config = dict(config) if config else {}
    config.setdefault("showLink", show_link)
    config.setdefault("linkText", link_text)

    # Get figure
    figure = tools.return_figure_from_figure_or_data(figure_or_data, validate)

    # Handle image request
    post_script = build_save_image_post_script(
        image, filename, image_height, image_width, "iplot"
    )

    # Show figure
    pio.show(
        figure,
        validate=validate,
        config=config,
        auto_play=auto_play,
        post_script=post_script,
        animation_opts=animation_opts,
    )


def plot(
    figure_or_data,
    show_link=False,
    link_text="Export to plot.ly",
    validate=True,
    output_type="file",
    include_plotlyjs=True,
    filename="temp-plot.html",
    auto_open=True,
    image=None,
    image_filename="plot_image",
    image_width=800,
    image_height=600,
    config=None,
    include_mathjax=False,
    auto_play=True,
    animation_opts=None,
):
    """Create a plotly graph locally as an HTML document or string.

    Example:
    ```
    from plotly.offline import plot
    import plotly.graph_objs as go

    plot([go.Scatter(x=[1, 2, 3], y=[3, 2, 6])], filename='my-graph.html')
    # We can also download an image of the plot by setting the image parameter
    # to the image format we want
    plot([go.Scatter(x=[1, 2, 3], y=[3, 2, 6])], filename='my-graph.html',
         image='jpeg')
    ```
    More examples below.

    figure_or_data -- a plotly.graph_objs.Figure or plotly.graph_objs.Data or
                      dict or list that describes a Plotly graph.
                      See https://plot.ly/python/ for examples of
                      graph descriptions.

    Keyword arguments:
    show_link (default=False) -- display a link in the bottom-right corner of
        of the chart that will export the chart to Plotly Cloud or
        Plotly Enterprise
    link_text (default='Export to plot.ly') -- the text of export link
    validate (default=True) -- validate that all of the keys in the figure
        are valid? omit if your version of plotly.js has become outdated
        with your version of graph_reference.json or if you need to include
        extra, unnecessary keys in your figure.
    output_type ('file' | 'div' - default 'file') -- if 'file', then
        the graph is saved as a standalone HTML file and `plot`
        returns None.
        If 'div', then `plot` returns a string that just contains the
        HTML <div> that contains the graph and the script to generate the
        graph.
        Use 'file' if you want to save and view a single graph at a time
        in a standalone HTML file.
        Use 'div' if you are embedding these graphs in an HTML file with
        other graphs or HTML markup, like a HTML report or an website.
    include_plotlyjs (True | False | 'cdn' | 'directory' | path - default=True)
        Specifies how the plotly.js library is included in the output html
        file or div string.

        If True, a script tag containing the plotly.js source code (~3MB)
        is included in the output.  HTML files generated with this option are
        fully self-contained and can be used offline.

        If 'cdn', a script tag that references the plotly.js CDN is included
        in the output. HTML files generated with this option are about 3MB
        smaller than those generated with include_plotlyjs=True, but they
        require an active internet connection in order to load the plotly.js
        library.

        If 'directory', a script tag is included that references an external
        plotly.min.js bundle that is assumed to reside in the same
        directory as the HTML file.  If output_type='file' then the
        plotly.min.js bundle is copied into the directory of the resulting
        HTML file. If a file named plotly.min.js already exists in the output
        directory then this file is left unmodified and no copy is performed.
        HTML files generated with this option can be used offline, but they
        require a copy of the plotly.min.js bundle in the same directory.
        This option is useful when many figures will be saved as HTML files in
        the same directory because the plotly.js source code will be included
        only once per output directory, rather than once per output file.

        If a string that ends in '.js', a script tag is included that
        references the specified path. This approach can be used to point
        the resulting HTML file to an alternative CDN.

        If False, no script tag referencing plotly.js is included. This is
        useful when output_type='div' and the resulting div string will be
        placed inside an HTML document that already loads plotly.js.  This
        option is not advised when output_type='file' as it will result in
        a non-functional html file.
    filename (default='temp-plot.html') -- The local filename to save the
        outputted chart to. If the filename already exists, it will be
        overwritten. This argument only applies if `output_type` is 'file'.
    auto_open (default=True) -- If True, open the saved file in a
        web browser after saving.
        This argument only applies if `output_type` is 'file'.
    image (default=None |'png' |'jpeg' |'svg' |'webp') -- This parameter sets
        the format of the image to be downloaded, if we choose to download an
        image. This parameter has a default value of None indicating that no
        image should be downloaded. Please note: for higher resolution images
        and more export options, consider making requests to our image servers.
        Type: `help(py.image)` for more details.
    image_filename (default='plot_image') -- Sets the name of the file your
        image will be saved to. The extension should not be included.
    image_height (default=600) -- Specifies the height of the image in `px`.
    image_width (default=800) -- Specifies the width of the image in `px`.
    config (default=None) -- Plot view options dictionary. Keyword arguments
        `show_link` and `link_text` set the associated options in this
        dictionary if it doesn't contain them already.
    include_mathjax (False | 'cdn' | path - default=False) --
        Specifies how the MathJax.js library is included in the output html
        file or div string.  MathJax is required in order to display labels
        with LaTeX typesetting.

        If False, no script tag referencing MathJax.js will be included in the
        output. HTML files generated with this option will not be able to
        display LaTeX typesetting.

        If 'cdn', a script tag that references a MathJax CDN location will be
        included in the output.  HTML files generated with this option will be
        able to display LaTeX typesetting as long as they have internet access.

        If a string that ends in '.js', a script tag is included that
        references the specified path. This approach can be used to point the
        resulting HTML file to an alternative CDN.
    auto_play (default=True) -- Whether to automatically start the animation
        sequence on page load if the figure contains frames. Has no effect if
        the figure does not contain frames.
    animation_opts (default=None) -- Dict of custom animation parameters that
        are used for the automatically started animation on page load. This
        dict is passed to the function Plotly.animate in Plotly.js. See
        https://github.com/plotly/plotly.js/blob/master/src/plots/animation_attributes.js
        for available options. Has no effect if the figure
        does not contain frames, or auto_play is False.

    Example:
    ```
    from plotly.offline import plot
    figure = {'data': [{'x': [0, 1], 'y': [0, 1]}],
              'layout': {'xaxis': {'range': [0, 5], 'autorange': False},
                         'yaxis': {'range': [0, 5], 'autorange': False},
                         'title': 'Start Title'},
              'frames': [{'data': [{'x': [1, 2], 'y': [1, 2]}]},
                         {'data': [{'x': [1, 4], 'y': [1, 4]}]},
                         {'data': [{'x': [3, 4], 'y': [3, 4]}],
                          'layout': {'title': 'End Title'}}]}
    plot(figure, animation_opts={'frame': {'duration': 1}})
    ```
    """
    import plotly.io as pio

    # Output type
    if output_type not in ["div", "file"]:
        raise ValueError(
            "`output_type` argument must be 'div' or 'file'. "
            "You supplied `" + output_type + "``"
        )
    if not filename.endswith(".html") and output_type == "file":
        warnings.warn(
            "Your filename `" + filename + "` didn't end with .html. "
            "Adding .html to the end of your file."
        )
        filename += ".html"

    # Config
    config = dict(config) if config else {}
    config.setdefault("showLink", show_link)
    config.setdefault("linkText", link_text)

    figure = tools.return_figure_from_figure_or_data(figure_or_data, validate)
    width = figure.get("layout", {}).get("width", "100%")
    height = figure.get("layout", {}).get("height", "100%")

    if width == "100%" or height == "100%":
        config.setdefault("responsive", True)

    # Handle image request
    post_script = build_save_image_post_script(
        image, image_filename, image_height, image_width, "plot"
    )

    if output_type == "file":
        pio.write_html(
            figure,
            filename,
            config=config,
            auto_play=auto_play,
            include_plotlyjs=include_plotlyjs,
            include_mathjax=include_mathjax,
            post_script=post_script,
            full_html=True,
            validate=validate,
            animation_opts=animation_opts,
            auto_open=auto_open,
        )
        return filename
    else:
        return pio.to_html(
            figure,
            config=config,
            auto_play=auto_play,
            include_plotlyjs=include_plotlyjs,
            include_mathjax=include_mathjax,
            post_script=post_script,
            full_html=False,
            validate=validate,
            animation_opts=animation_opts,
        )


def plot_mpl(
    mpl_fig,
    resize=False,
    strip_style=False,
    verbose=False,
    show_link=False,
    link_text="Export to plot.ly",
    validate=True,
    output_type="file",
    include_plotlyjs=True,
    filename="temp-plot.html",
    auto_open=True,
    image=None,
    image_filename="plot_image",
    image_height=600,
    image_width=800,
):
    """
    Convert a matplotlib figure to a Plotly graph stored locally as HTML.

    For more information on converting matplotlib visualizations to plotly
    graphs, call help(plotly.tools.mpl_to_plotly)

    For more information on creating plotly charts locally as an HTML document
    or string, call help(plotly.offline.plot)

    mpl_fig -- a matplotlib figure object to convert to a plotly graph

    Keyword arguments:
    resize (default=False) -- allow plotly to choose the figure size.
    strip_style (default=False) -- allow plotly to choose style options.
    verbose (default=False) -- print message.
    show_link (default=False) -- display a link in the bottom-right corner of
        of the chart that will export the chart to Plotly Cloud or
        Plotly Enterprise
    link_text (default='Export to plot.ly') -- the text of export link
    validate (default=True) -- validate that all of the keys in the figure
        are valid? omit if your version of plotly.js has become outdated
        with your version of graph_reference.json or if you need to include
        extra, unnecessary keys in your figure.
    output_type ('file' | 'div' - default 'file') -- if 'file', then
        the graph is saved as a standalone HTML file and `plot`
        returns None.
        If 'div', then `plot` returns a string that just contains the
        HTML <div> that contains the graph and the script to generate the
        graph.
        Use 'file' if you want to save and view a single graph at a time
        in a standalone HTML file.
        Use 'div' if you are embedding these graphs in an HTML file with
        other graphs or HTML markup, like a HTML report or an website.
    include_plotlyjs (default=True) -- If True, include the plotly.js
        source code in the output file or string.
        Set as False if your HTML file already contains a copy of the plotly.js
        library.
    filename (default='temp-plot.html') -- The local filename to save the
        outputted chart to. If the filename already exists, it will be
        overwritten. This argument only applies if `output_type` is 'file'.
    auto_open (default=True) -- If True, open the saved file in a
        web browser after saving.
        This argument only applies if `output_type` is 'file'.
    image (default=None |'png' |'jpeg' |'svg' |'webp') -- This parameter sets
        the format of the image to be downloaded, if we choose to download an
        image. This parameter has a default value of None indicating that no
        image should be downloaded.
    image_filename (default='plot_image') -- Sets the name of the file your
        image will be saved to. The extension should not be included.
    image_height (default=600) -- Specifies the height of the image in `px`.
    image_width (default=800) -- Specifies the width of the image in `px`.

    Example:
    ```
    from plotly.offline import init_notebook_mode, plot_mpl
    import matplotlib.pyplot as plt

    init_notebook_mode()

    fig = plt.figure()
    x = [10, 15, 20, 25, 30]
    y = [100, 250, 200, 150, 300]
    plt.plot(x, y, "o")

    plot_mpl(fig)
    # If you want to to download an image of the figure as well
    plot_mpl(fig, image='png')
    ```
    """
    plotly_plot = tools.mpl_to_plotly(mpl_fig, resize, strip_style, verbose)
    return plot(
        plotly_plot,
        show_link,
        link_text,
        validate,
        output_type,
        include_plotlyjs,
        filename,
        auto_open,
        image=image,
        image_filename=image_filename,
        image_height=image_height,
        image_width=image_width,
    )


def iplot_mpl(
    mpl_fig,
    resize=False,
    strip_style=False,
    verbose=False,
    show_link=False,
    link_text="Export to plot.ly",
    validate=True,
    image=None,
    image_filename="plot_image",
    image_height=600,
    image_width=800,
):
    """
    Convert a matplotlib figure to a plotly graph and plot inside an IPython
    notebook without connecting to an external server.

    To save the chart to Plotly Cloud or Plotly Enterprise, use
    `plotly.plotly.plot_mpl`.

    For more information on converting matplotlib visualizations to plotly
    graphs call `help(plotly.tools.mpl_to_plotly)`

    For more information on plotting plotly charts offline in an Ipython
    notebook call `help(plotly.offline.iplot)`

    mpl_fig -- a matplotlib.figure to convert to a plotly graph

    Keyword arguments:
    resize (default=False) -- allow plotly to choose the figure size.
    strip_style (default=False) -- allow plotly to choose style options.
    verbose (default=False) -- print message.
    show_link (default=False) -- display a link in the bottom-right corner of
                                of the chart that will export the chart to
                                Plotly Cloud or Plotly Enterprise
    link_text (default='Export to plot.ly') -- the text of export link
    validate (default=True) -- validate that all of the keys in the figure
                               are valid? omit if your version of plotly.js
                               has become outdated with your version of
                               graph_reference.json or if you need to include
                               extra, unnecessary keys in your figure.
    image (default=None |'png' |'jpeg' |'svg' |'webp') -- This parameter sets
        the format of the image to be downloaded, if we choose to download an
        image. This parameter has a default value of None indicating that no
        image should be downloaded.
    image_filename (default='plot_image') -- Sets the name of the file your
        image will be saved to. The extension should not be included.
    image_height (default=600) -- Specifies the height of the image in `px`.
    image_width (default=800) -- Specifies the width of the image in `px`.

    Example:
    ```
    from plotly.offline import init_notebook_mode, iplot_mpl
    import matplotlib.pyplot as plt

    fig = plt.figure()
    x = [10, 15, 20, 25, 30]
    y = [100, 250, 200, 150, 300]
    plt.plot(x, y, "o")

    init_notebook_mode()
    iplot_mpl(fig)
    # and if you want to download an image of the figure as well
    iplot_mpl(fig, image='jpeg')
    ```
    """
    plotly_plot = tools.mpl_to_plotly(mpl_fig, resize, strip_style, verbose)
    return iplot(
        plotly_plot,
        show_link,
        link_text,
        validate,
        image=image,
        filename=image_filename,
        image_height=image_height,
        image_width=image_width,
    )


def enable_mpl_offline(
    resize=False,
    strip_style=False,
    verbose=False,
    show_link=False,
    link_text="Export to plot.ly",
    validate=True,
):
    """
    Convert mpl plots to locally hosted HTML documents.

    This function should be used with the inline matplotlib backend
    that ships with IPython that can be enabled with `%pylab inline`
    or `%matplotlib inline`. This works by adding an HTML formatter
    for Figure objects; the existing SVG/PNG formatters will remain
    enabled.

    (idea taken from `mpld3._display.enable_notebook`)

    Example:
    ```
    from plotly.offline import enable_mpl_offline
    import matplotlib.pyplot as plt

    enable_mpl_offline()

    fig = plt.figure()
    x = [10, 15, 20, 25, 30]
    y = [100, 250, 200, 150, 300]
    plt.plot(x, y, "o")
    fig
    ```
    """
    init_notebook_mode()
    ipython = get_module("IPython")
    matplotlib = get_module("matplotlib")

    ip = ipython.core.getipython.get_ipython()
    formatter = ip.display_formatter.formatters["text/html"]
    formatter.for_type(
        matplotlib.figure.Figure,
        lambda fig: iplot_mpl(
            fig, resize, strip_style, verbose, show_link, link_text, validate
        ),
    )
