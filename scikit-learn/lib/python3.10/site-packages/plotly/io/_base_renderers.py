import base64
import json
import webbrowser
import inspect
import os
from os.path import isdir

from plotly import optional_imports
from plotly.io import to_json, to_image, write_image, write_html
from plotly.io._utils import plotly_cdn_url
from plotly.offline.offline import _get_jconfig, get_plotlyjs
from plotly.tools import return_figure_from_figure_or_data

ipython_display = optional_imports.get_module("IPython.display")
IPython = optional_imports.get_module("IPython")

try:
    from http.server import BaseHTTPRequestHandler, HTTPServer
except ImportError:
    # Python 2.7
    from BaseHTTPServer import BaseHTTPRequestHandler, HTTPServer


class BaseRenderer(object):
    """
    Base class for all renderers
    """

    def activate(self):
        pass

    def __repr__(self):
        try:
            init_sig = inspect.signature(self.__init__)
            init_args = list(init_sig.parameters.keys())
        except AttributeError:
            # Python 2.7
            argspec = inspect.getargspec(self.__init__)
            init_args = [a for a in argspec.args if a != "self"]

        return "{cls}({attrs})\n{doc}".format(
            cls=self.__class__.__name__,
            attrs=", ".join("{}={!r}".format(k, self.__dict__[k]) for k in init_args),
            doc=self.__doc__,
        )

    def __hash__(self):
        # Constructor args fully define uniqueness
        return hash(repr(self))


class MimetypeRenderer(BaseRenderer):
    """
    Base class for all mime type renderers
    """

    def to_mimebundle(self, fig_dict):
        raise NotImplementedError()


class JsonRenderer(MimetypeRenderer):
    """
    Renderer to display figures as JSON hierarchies.  This renderer is
    compatible with JupyterLab and VSCode.

    mime type: 'application/json'
    """

    def to_mimebundle(self, fig_dict):
        value = json.loads(to_json(fig_dict, validate=False, remove_uids=False))
        return {"application/json": value}


# Plotly mimetype
class PlotlyRenderer(MimetypeRenderer):
    """
    Renderer to display figures using the plotly mime type.  This renderer is
    compatible with VSCode and nteract.

    mime type: 'application/vnd.plotly.v1+json'
    """

    def __init__(self, config=None):
        self.config = dict(config) if config else {}

    def to_mimebundle(self, fig_dict):
        config = _get_jconfig(self.config)
        if config:
            fig_dict["config"] = config

        json_compatible_fig_dict = json.loads(
            to_json(fig_dict, validate=False, remove_uids=False)
        )

        return {"application/vnd.plotly.v1+json": json_compatible_fig_dict}


# Static Image
class ImageRenderer(MimetypeRenderer):
    """
    Base class for all static image renderers
    """

    def __init__(
        self,
        mime_type,
        b64_encode=False,
        format=None,
        width=None,
        height=None,
        scale=None,
        engine="auto",
    ):
        self.mime_type = mime_type
        self.b64_encode = b64_encode
        self.format = format
        self.width = width
        self.height = height
        self.scale = scale
        self.engine = engine

    def to_mimebundle(self, fig_dict):
        image_bytes = to_image(
            fig_dict,
            format=self.format,
            width=self.width,
            height=self.height,
            scale=self.scale,
            validate=False,
            engine=self.engine,
        )

        if self.b64_encode:
            image_str = base64.b64encode(image_bytes).decode("utf8")
        else:
            image_str = image_bytes.decode("utf8")

        return {self.mime_type: image_str}


class PngRenderer(ImageRenderer):
    """
    Renderer to display figures as static PNG images.  This renderer requires
    either the kaleido package or the orca command-line utility and is broadly
    compatible across IPython environments (classic Jupyter Notebook, JupyterLab,
    QtConsole, VSCode, PyCharm, etc) and nbconvert targets (HTML, PDF, etc.).

    mime type: 'image/png'
    """

    def __init__(self, width=None, height=None, scale=None, engine=None):
        super(PngRenderer, self).__init__(
            mime_type="image/png",
            b64_encode=True,
            format="png",
            width=width,
            height=height,
            scale=scale,
            engine=engine,
        )


class SvgRenderer(ImageRenderer):
    """
    Renderer to display figures as static SVG images.  This renderer requires
    either the kaleido package or the orca command-line utility and is broadly
    compatible across IPython environments (classic Jupyter Notebook, JupyterLab,
    QtConsole, VSCode, PyCharm, etc) and nbconvert targets (HTML, PDF, etc.).

    mime type: 'image/svg+xml'
    """

    def __init__(self, width=None, height=None, scale=None, engine=None):
        super(SvgRenderer, self).__init__(
            mime_type="image/svg+xml",
            b64_encode=False,
            format="svg",
            width=width,
            height=height,
            scale=scale,
            engine=engine,
        )


class JpegRenderer(ImageRenderer):
    """
    Renderer to display figures as static JPEG images.  This renderer requires
    either the kaleido package or the orca command-line utility and is broadly
    compatible across IPython environments (classic Jupyter Notebook, JupyterLab,
    QtConsole, VSCode, PyCharm, etc) and nbconvert targets (HTML, PDF, etc.).

    mime type: 'image/jpeg'
    """

    def __init__(self, width=None, height=None, scale=None, engine=None):
        super(JpegRenderer, self).__init__(
            mime_type="image/jpeg",
            b64_encode=True,
            format="jpg",
            width=width,
            height=height,
            scale=scale,
            engine=engine,
        )


class PdfRenderer(ImageRenderer):
    """
    Renderer to display figures as static PDF images.  This renderer requires
    either the kaleido package or the orca command-line utility and is compatible
    with JupyterLab and the LaTeX-based nbconvert export to PDF.

    mime type: 'application/pdf'
    """

    def __init__(self, width=None, height=None, scale=None, engine=None):
        super(PdfRenderer, self).__init__(
            mime_type="application/pdf",
            b64_encode=True,
            format="pdf",
            width=width,
            height=height,
            scale=scale,
            engine=engine,
        )


# HTML
# Build script to set global PlotlyConfig object. This must execute before
# plotly.js is loaded.
_window_plotly_config = """\
window.PlotlyConfig = {MathJaxConfig: 'local'};"""

_mathjax_config = """\
if (window.MathJax && window.MathJax.Hub && window.MathJax.Hub.Config) {window.MathJax.Hub.Config({SVG: {font: "STIX-Web"}});}"""


class HtmlRenderer(MimetypeRenderer):
    """
    Base class for all HTML mime type renderers

    mime type: 'text/html'
    """

    def __init__(
        self,
        connected=False,
        full_html=False,
        global_init=False,
        config=None,
        auto_play=False,
        post_script=None,
        animation_opts=None,
        include_plotlyjs=True,
    ):
        self.config = dict(config) if config else {}
        self.auto_play = auto_play
        self.connected = connected
        self.global_init = global_init
        self.full_html = full_html
        self.animation_opts = animation_opts
        self.post_script = post_script
        self.include_plotlyjs = "cdn" if self.connected else include_plotlyjs

    def activate(self):
        if self.global_init:
            if not ipython_display:
                raise ValueError(
                    "The {cls} class requires ipython but it is not installed".format(
                        cls=self.__class__.__name__
                    )
                )

            if self.connected:
                script = """\
        <script type="text/javascript">
        {win_config}
        {mathjax_config}
        </script>
        <script type="module">import \"{plotly_cdn}\"</script>
        """.format(
                    win_config=_window_plotly_config,
                    mathjax_config=_mathjax_config,
                    plotly_cdn=plotly_cdn_url().rstrip(".js"),
                )

            else:
                # If not connected then we embed a copy of the plotly.js
                # library in the notebook
                script = """\
        <script type="text/javascript">
        {win_config}
        {mathjax_config}
        </script>
        <script>{script}</script>
        """.format(
                    script=get_plotlyjs(),
                    win_config=_window_plotly_config,
                    mathjax_config=_mathjax_config,
                )

            ipython_display.display_html(script, raw=True)

    def to_mimebundle(self, fig_dict):
        from plotly.io import to_html

        include_mathjax = "cdn"

        # build post script
        post_script = [
            """
var gd = document.getElementById('{plot_id}');
var x = new MutationObserver(function (mutations, observer) {{
        var display = window.getComputedStyle(gd).display;
        if (!display || display === 'none') {{
            console.log([gd, 'removed!']);
            Plotly.purge(gd);
            observer.disconnect();
        }}
}});

// Listen for the removal of the full notebook cells
var notebookContainer = gd.closest('#notebook-container');
if (notebookContainer) {{
    x.observe(notebookContainer, {childList: true});
}}

// Listen for the clearing of the current output cell
var outputEl = gd.closest('.output');
if (outputEl) {{
    x.observe(outputEl, {childList: true});
}}
"""
        ]

        # Add user defined post script
        if self.post_script:
            if not isinstance(self.post_script, (list, tuple)):
                post_script.append(self.post_script)
            else:
                post_script.extend(self.post_script)

        html = to_html(
            fig_dict,
            config=self.config,
            auto_play=self.auto_play,
            include_plotlyjs=self.include_plotlyjs,
            include_mathjax=include_mathjax,
            post_script=post_script,
            full_html=self.full_html,
            animation_opts=self.animation_opts,
            default_width="100%",
            default_height=525,
            validate=False,
        )

        return {"text/html": html}


class NotebookRenderer(HtmlRenderer):
    """
    Renderer to display interactive figures in the classic Jupyter Notebook.
    This renderer is also useful for notebooks that will be converted to
    HTML using nbconvert/nbviewer as it will produce standalone HTML files
    that include interactive figures.

    This renderer automatically performs global notebook initialization when
    activated.

    mime type: 'text/html'
    """

    def __init__(
        self,
        connected=False,
        config=None,
        auto_play=False,
        post_script=None,
        animation_opts=None,
        include_plotlyjs=False,
    ):
        super(NotebookRenderer, self).__init__(
            connected=connected,
            full_html=False,
            global_init=True,
            config=config,
            auto_play=auto_play,
            post_script=post_script,
            animation_opts=animation_opts,
            include_plotlyjs=include_plotlyjs,
        )


class KaggleRenderer(HtmlRenderer):
    """
    Renderer to display interactive figures in Kaggle Notebooks.

    Same as NotebookRenderer but with connected=True so that the plotly.js
    bundle is loaded from a CDN rather than being embedded in the notebook.

    This renderer is enabled by default when running in a Kaggle notebook.

    mime type: 'text/html'
    """

    def __init__(
        self, config=None, auto_play=False, post_script=None, animation_opts=None
    ):
        super(KaggleRenderer, self).__init__(
            connected=True,
            full_html=False,
            global_init=True,
            config=config,
            auto_play=auto_play,
            post_script=post_script,
            animation_opts=animation_opts,
            include_plotlyjs=False,
        )


class AzureRenderer(HtmlRenderer):
    """
    Renderer to display interactive figures in Azure Notebooks.

    Same as NotebookRenderer but with connected=True so that the plotly.js
    bundle is loaded from a CDN rather than being embedded in the notebook.

    This renderer is enabled by default when running in an Azure notebook.

    mime type: 'text/html'
    """

    def __init__(
        self, config=None, auto_play=False, post_script=None, animation_opts=None
    ):
        super(AzureRenderer, self).__init__(
            connected=True,
            full_html=False,
            global_init=True,
            config=config,
            auto_play=auto_play,
            post_script=post_script,
            animation_opts=animation_opts,
            include_plotlyjs=False,
        )


class ColabRenderer(HtmlRenderer):
    """
    Renderer to display interactive figures in Google Colab Notebooks.

    This renderer is enabled by default when running in a Colab notebook.

    mime type: 'text/html'
    """

    def __init__(
        self, config=None, auto_play=False, post_script=None, animation_opts=None
    ):
        super(ColabRenderer, self).__init__(
            connected=True,
            full_html=True,
            global_init=False,
            config=config,
            auto_play=auto_play,
            post_script=post_script,
            animation_opts=animation_opts,
        )


class IFrameRenderer(MimetypeRenderer):
    """
    Renderer to display interactive figures using an IFrame.  HTML
    representations of Figures are saved to an `iframe_figures/` directory and
    iframe HTML elements that reference these files are inserted into the
    notebook.

    With this approach, neither plotly.js nor the figure data are embedded in
    the notebook, so this is a good choice for notebooks that contain so many
    large figures that basic operations (like saving and opening) become
    very slow.

    Notebooks using this renderer will display properly when exported to HTML
    as long as the `iframe_figures/` directory is placed in the same directory
    as the exported html file.

    Note that the HTML files in `iframe_figures/` are numbered according to
    the IPython cell execution count and so they will start being overwritten
    each time the kernel is restarted.  This directory may be deleted whenever
    the kernel is restarted and it will be automatically recreated.

    mime type: 'text/html'
    """

    def __init__(
        self,
        config=None,
        auto_play=False,
        post_script=None,
        animation_opts=None,
        include_plotlyjs=True,
        html_directory="iframe_figures",
    ):
        self.config = config
        self.auto_play = auto_play
        self.post_script = post_script
        self.animation_opts = animation_opts
        self.include_plotlyjs = include_plotlyjs
        self.html_directory = html_directory

    def to_mimebundle(self, fig_dict):
        from plotly.io import write_html

        # Make iframe size slightly larger than figure size to avoid
        # having iframe have its own scroll bar.
        iframe_buffer = 20
        layout = fig_dict.get("layout", {})

        if layout.get("width", False):
            iframe_width = str(layout["width"] + iframe_buffer) + "px"
        else:
            iframe_width = "100%"

        if layout.get("height", False):
            iframe_height = layout["height"] + iframe_buffer
        else:
            iframe_height = str(525 + iframe_buffer) + "px"

        # Build filename using ipython cell number
        filename = self.build_filename()

        # Make directory for
        try:
            os.makedirs(self.html_directory)
        except OSError:
            if not isdir(self.html_directory):
                raise

        write_html(
            fig_dict,
            filename,
            config=self.config,
            auto_play=self.auto_play,
            include_plotlyjs=self.include_plotlyjs,
            include_mathjax="cdn",
            auto_open=False,
            post_script=self.post_script,
            animation_opts=self.animation_opts,
            default_width="100%",
            default_height=525,
            validate=False,
        )

        # Build IFrame
        iframe_html = """\
<iframe
    scrolling="no"
    width="{width}"
    height="{height}"
    src="{src}"
    frameborder="0"
    allowfullscreen
></iframe>
""".format(width=iframe_width, height=iframe_height, src=self.build_url(filename))

        return {"text/html": iframe_html}

    def build_filename(self):
        ip = IPython.get_ipython() if IPython else None
        try:
            cell_number = list(ip.history_manager.get_tail(1))[0][1] + 1 if ip else 0
        except Exception:
            cell_number = 0
        return "{dirname}/figure_{cell_number}.html".format(
            dirname=self.html_directory, cell_number=cell_number
        )

    def build_url(self, filename):
        return filename


class CoCalcRenderer(IFrameRenderer):
    _render_count = 0

    def build_filename(self):
        filename = "{dirname}/figure_{render_count}.html".format(
            dirname=self.html_directory, render_count=CoCalcRenderer._render_count
        )

        CoCalcRenderer._render_count += 1
        return filename

    def build_url(self, filename):
        return "{filename}?fullscreen=kiosk".format(filename=filename)


class ExternalRenderer(BaseRenderer):
    """
    Base class for external renderers.  ExternalRenderer subclasses
    do not display figures inline in a notebook environment, but render
    figures by some external means (e.g. a separate browser tab).

    Unlike MimetypeRenderer subclasses, ExternalRenderer subclasses are not
    invoked when a figure is asked to display itself in the notebook.
    Instead, they are invoked when the plotly.io.show function is called
    on a figure.
    """

    def render(self, fig):
        raise NotImplementedError()


def open_html_in_browser(html, using=None, new=0, autoraise=True):
    """
    Display html in a web browser without creating a temp file.

    Instantiates a trivial http server and uses the webbrowser module to
    open a URL to retrieve html from that server.

    Parameters
    ----------
    html: str
        HTML string to display
    using, new, autoraise:
        See docstrings in webbrowser.get and webbrowser.open
    """
    if isinstance(html, str):
        html = html.encode("utf8")

    browser = None

    if using is None:
        browser = webbrowser.get(None)
    else:
        if not isinstance(using, tuple):
            using = (using,)
        for browser_key in using:
            try:
                browser = webbrowser.get(browser_key)
                if browser is not None:
                    break
            except webbrowser.Error:
                pass

        if browser is None:
            raise ValueError("Can't locate a browser with key in " + str(using))

    class OneShotRequestHandler(BaseHTTPRequestHandler):
        def do_GET(self):
            self.send_response(200)
            self.send_header("Content-type", "text/html")
            self.end_headers()

            bufferSize = 1024 * 1024
            for i in range(0, len(html), bufferSize):
                self.wfile.write(html[i : i + bufferSize])

        def log_message(self, format, *args):
            # Silence stderr logging
            pass

    server = HTTPServer(("127.0.0.1", 0), OneShotRequestHandler)
    browser.open(
        "http://127.0.0.1:%s" % server.server_port, new=new, autoraise=autoraise
    )

    server.handle_request()


class BrowserRenderer(ExternalRenderer):
    """
    Renderer to display interactive figures in an external web browser.
    This renderer will open a new browser window or tab when the
    plotly.io.show function is called on a figure.

    This renderer has no ipython/jupyter dependencies and is a good choice
    for use in environments that do not support the inline display of
    interactive figures.

    mime type: 'text/html'
    """

    def __init__(
        self,
        config=None,
        auto_play=False,
        using=None,
        new=0,
        autoraise=True,
        post_script=None,
        animation_opts=None,
    ):
        self.config = config
        self.auto_play = auto_play
        self.using = using
        self.new = new
        self.autoraise = autoraise
        self.post_script = post_script
        self.animation_opts = animation_opts

    def render(self, fig_dict):
        from plotly.io import to_html

        html = to_html(
            fig_dict,
            config=self.config,
            auto_play=self.auto_play,
            include_plotlyjs=True,
            include_mathjax="cdn",
            post_script=self.post_script,
            full_html=True,
            animation_opts=self.animation_opts,
            default_width="100%",
            default_height="100%",
            validate=False,
        )
        open_html_in_browser(html, self.using, self.new, self.autoraise)


class DatabricksRenderer(ExternalRenderer):
    def __init__(
        self,
        config=None,
        auto_play=False,
        post_script=None,
        animation_opts=None,
        include_plotlyjs="cdn",
    ):
        self.config = config
        self.auto_play = auto_play
        self.post_script = post_script
        self.animation_opts = animation_opts
        self.include_plotlyjs = include_plotlyjs
        self._displayHTML = None

    @property
    def displayHTML(self):
        import inspect

        if self._displayHTML is None:
            for frame in inspect.getouterframes(inspect.currentframe()):
                global_names = set(frame.frame.f_globals)
                # Check for displayHTML plus a few others to reduce chance of a false
                # hit.
                if all(v in global_names for v in ["displayHTML", "display", "spark"]):
                    self._displayHTML = frame.frame.f_globals["displayHTML"]
                    break

            if self._displayHTML is None:
                raise EnvironmentError(
                    """
Unable to detect the Databricks displayHTML function. The 'databricks' renderer is only
supported when called from within the Databricks notebook environment."""
                )

        return self._displayHTML

    def render(self, fig_dict):
        from plotly.io import to_html

        html = to_html(
            fig_dict,
            config=self.config,
            auto_play=self.auto_play,
            include_plotlyjs=self.include_plotlyjs,
            include_mathjax="cdn",
            post_script=self.post_script,
            full_html=True,
            animation_opts=self.animation_opts,
            default_width="100%",
            default_height="100%",
            validate=False,
        )

        # displayHTML is a Databricks notebook built-in function
        self.displayHTML(html)


class SphinxGalleryHtmlRenderer(HtmlRenderer):
    def __init__(
        self,
        connected=True,
        config=None,
        auto_play=False,
        post_script=None,
        animation_opts=None,
    ):
        super(SphinxGalleryHtmlRenderer, self).__init__(
            connected=connected,
            full_html=False,
            global_init=False,
            config=config,
            auto_play=auto_play,
            post_script=post_script,
            animation_opts=animation_opts,
        )

    def to_mimebundle(self, fig_dict):
        from plotly.io import to_html

        if self.connected:
            include_plotlyjs = "cdn"
            include_mathjax = "cdn"
        else:
            include_plotlyjs = True
            include_mathjax = "cdn"

        html = to_html(
            fig_dict,
            config=self.config,
            auto_play=self.auto_play,
            include_plotlyjs=include_plotlyjs,
            include_mathjax=include_mathjax,
            full_html=self.full_html,
            animation_opts=self.animation_opts,
            default_width="100%",
            default_height=525,
            validate=False,
        )

        return {"text/html": html}


class SphinxGalleryOrcaRenderer(ExternalRenderer):
    def render(self, fig_dict):
        stack = inspect.stack()
        # Name of script from which plot function was called is retrieved
        try:
            filename = stack[3].filename  # let's hope this is robust...
        except Exception:  # python 2
            filename = stack[3][1]
        filename_root, _ = os.path.splitext(filename)
        filename_html = filename_root + ".html"
        filename_png = filename_root + ".png"
        figure = return_figure_from_figure_or_data(fig_dict, True)
        _ = write_html(fig_dict, file=filename_html, include_plotlyjs="cdn")
        try:
            write_image(figure, filename_png)
        except (ValueError, ImportError):
            raise ImportError(
                "orca and psutil are required to use the `sphinx-gallery-orca` renderer. "
                "See https://plotly.com/python/static-image-export/ for instructions on "
                "how to install orca. Alternatively, you can use the `sphinx-gallery` "
                "renderer (note that png thumbnails can only be generated with "
                "the `sphinx-gallery-orca` renderer)."
            )
