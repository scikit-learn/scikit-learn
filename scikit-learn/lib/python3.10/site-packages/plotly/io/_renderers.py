import textwrap
from copy import copy
import os
from packaging.version import Version
import warnings

from plotly import optional_imports

from plotly.io._base_renderers import (
    MimetypeRenderer,
    ExternalRenderer,
    PlotlyRenderer,
    NotebookRenderer,
    KaggleRenderer,
    AzureRenderer,
    ColabRenderer,
    JsonRenderer,
    PngRenderer,
    JpegRenderer,
    SvgRenderer,
    PdfRenderer,
    BrowserRenderer,
    IFrameRenderer,
    SphinxGalleryHtmlRenderer,
    SphinxGalleryOrcaRenderer,
    CoCalcRenderer,
    DatabricksRenderer,
)
from plotly.io._utils import validate_coerce_fig_to_dict

ipython = optional_imports.get_module("IPython")
ipython_display = optional_imports.get_module("IPython.display")
nbformat = optional_imports.get_module("nbformat")


def display_jupyter_version_warnings():
    parent_process = None
    try:
        psutil = optional_imports.get_module("psutil")
        if psutil is not None:
            parent_process = psutil.Process().parent().cmdline()[-1]
    except Exception:
        pass

    if parent_process is None:
        return
    elif "jupyter-notebook" in parent_process:
        jupyter_notebook = optional_imports.get_module("notebook")
        if jupyter_notebook is not None and jupyter_notebook.__version__ < "7":
            # Add warning about upgrading notebook
            warnings.warn(
                f"Plotly version >= 6 requires Jupyter Notebook >= 7 but you have {jupyter_notebook.__version__} installed.\n To upgrade Jupyter Notebook, please run `pip install notebook --upgrade`."
            )
    elif "jupyter-lab" in parent_process:
        jupyter_lab = optional_imports.get_module("jupyterlab")
        if jupyter_lab is not None and jupyter_lab.__version__ < "3":
            # Add warning about upgrading jupyterlab
            warnings.warn(
                f"Plotly version >= 6 requires JupyterLab >= 3 but you have {jupyter_lab.__version__} installed. To upgrade JupyterLab, please run `pip install jupyterlab --upgrade`."
            )


# Renderer configuration class
# -----------------------------
class RenderersConfig(object):
    """
    Singleton object containing the current renderer configurations
    """

    def __init__(self):
        self._renderers = {}
        self._default_name = None
        self._default_renderers = []
        self._render_on_display = False
        self._to_activate = []

    # ### Magic methods ###
    # Make this act as a dict of renderers
    def __len__(self):
        return len(self._renderers)

    def __contains__(self, item):
        return item in self._renderers

    def __iter__(self):
        return iter(self._renderers)

    def __getitem__(self, item):
        renderer = self._renderers[item]
        return renderer

    def __setitem__(self, key, value):
        if not isinstance(value, (MimetypeRenderer, ExternalRenderer)):
            raise ValueError(
                """\
Renderer must be a subclass of MimetypeRenderer or ExternalRenderer.
    Received value with type: {typ}""".format(typ=type(value))
            )

        self._renderers[key] = value

    def __delitem__(self, key):
        # Remove template
        del self._renderers[key]

        # Check if we need to remove it as the default
        if self._default == key:
            self._default = None

    def keys(self):
        return self._renderers.keys()

    def items(self):
        return self._renderers.items()

    def update(self, d={}, **kwargs):
        """
        Update one or more renderers from a dict or from input keyword
        arguments.

        Parameters
        ----------
        d: dict
            Dictionary from renderer names to new renderer objects.

        kwargs
            Named argument value pairs where the name is a renderer name
            and the value is a new renderer object
        """
        for k, v in dict(d, **kwargs).items():
            self[k] = v

    # ### Properties ###
    @property
    def default(self):
        """
        The default renderer, or None if no there is no default

        If not None, the default renderer is used to render
        figures when the `plotly.io.show` function is called on a Figure.

        If `plotly.io.renderers.render_on_display` is True, then the default
        renderer will also be used to display Figures automatically when
        displayed in the Jupyter Notebook

        Multiple renderers may be registered by separating their names with
        '+' characters. For example, to specify rendering compatible with
        the classic Jupyter Notebook, JupyterLab, and PDF export:

        >>> import plotly.io as pio
        >>> pio.renderers.default = 'notebook+jupyterlab+pdf'

        The names of available renderers may be retrieved with:

        >>> import plotly.io as pio
        >>> list(pio.renderers)

        Returns
        -------
        str
        """
        return self._default_name

    @default.setter
    def default(self, value):
        # Handle None
        if not value:
            # _default_name should always be a string so we can do
            # pio.renderers.default.split('+')
            self._default_name = ""
            self._default_renderers = []
            return

        # Store defaults name and list of renderer(s)
        renderer_names = self._validate_coerce_renderers(value)
        self._default_name = value
        self._default_renderers = [self[name] for name in renderer_names]

        # Register renderers for activation before their next use
        self._to_activate = list(self._default_renderers)

    @property
    def render_on_display(self):
        """
        If True, the default mimetype renderers will be used to render
        figures when they are displayed in an IPython context.

        Returns
        -------
        bool
        """
        return self._render_on_display

    @render_on_display.setter
    def render_on_display(self, val):
        self._render_on_display = bool(val)

    def _activate_pending_renderers(self, cls=object):
        """
        Activate all renderers that are waiting in the _to_activate list

        Parameters
        ----------
        cls
            Only activate renders that are subclasses of this class
        """
        to_activate_with_cls = [
            r for r in self._to_activate if cls and isinstance(r, cls)
        ]

        while to_activate_with_cls:
            # Activate renderers from left to right so that right-most
            # renderers take precedence
            renderer = to_activate_with_cls.pop(0)
            renderer.activate()

        self._to_activate = [
            r for r in self._to_activate if not (cls and isinstance(r, cls))
        ]

    def _validate_coerce_renderers(self, renderers_string):
        """
        Input a string and validate that it contains the names of one or more
        valid renderers separated on '+' characters.  If valid, return
        a list of the renderer names

        Parameters
        ----------
        renderers_string: str

        Returns
        -------
        list of str
        """
        # Validate value
        if not isinstance(renderers_string, str):
            raise ValueError("Renderer must be specified as a string")

        renderer_names = renderers_string.split("+")
        invalid = [name for name in renderer_names if name not in self]
        if invalid:
            raise ValueError(
                """
Invalid named renderer(s) received: {}""".format(str(invalid))
            )

        return renderer_names

    def __repr__(self):
        return """\
Renderers configuration
-----------------------
    Default renderer: {default}
    Available renderers:
{available}
""".format(default=repr(self.default), available=self._available_renderers_str())

    def _available_renderers_str(self):
        """
        Return nicely wrapped string representation of all
        available renderer names
        """
        available = "\n".join(
            textwrap.wrap(
                repr(list(self)),
                width=79 - 8,
                initial_indent=" " * 8,
                subsequent_indent=" " * 9,
            )
        )
        return available

    def _build_mime_bundle(self, fig_dict, renderers_string=None, **kwargs):
        """
        Build a mime bundle dict containing a kev/value pair for each
        MimetypeRenderer specified in either the default renderer string,
        or in the supplied renderers_string argument.

        Note that this method skips any renderers that are not subclasses
        of MimetypeRenderer.

        Parameters
        ----------
        fig_dict: dict
            Figure dictionary
        renderers_string: str or None (default None)
            Renderer string to process rather than the current default
            renderer string

        Returns
        -------
        dict
        """
        if renderers_string:
            renderer_names = self._validate_coerce_renderers(renderers_string)
            renderers_list = [self[name] for name in renderer_names]

            # Activate these non-default renderers
            for renderer in renderers_list:
                if isinstance(renderer, MimetypeRenderer):
                    renderer.activate()
        else:
            # Activate any pending default renderers
            self._activate_pending_renderers(cls=MimetypeRenderer)
            renderers_list = self._default_renderers

        bundle = {}
        for renderer in renderers_list:
            if isinstance(renderer, MimetypeRenderer):
                renderer = copy(renderer)
                for k, v in kwargs.items():
                    if hasattr(renderer, k):
                        setattr(renderer, k, v)

                bundle.update(renderer.to_mimebundle(fig_dict))

        return bundle

    def _perform_external_rendering(self, fig_dict, renderers_string=None, **kwargs):
        """
        Perform external rendering for each ExternalRenderer specified
        in either the default renderer string, or in the supplied
        renderers_string argument.

        Note that this method skips any renderers that are not subclasses
        of ExternalRenderer.

        Parameters
        ----------
        fig_dict: dict
            Figure dictionary
        renderers_string: str or None (default None)
            Renderer string to process rather than the current default
            renderer string

        Returns
        -------
        None
        """
        if renderers_string:
            renderer_names = self._validate_coerce_renderers(renderers_string)
            renderers_list = [self[name] for name in renderer_names]

            # Activate these non-default renderers
            for renderer in renderers_list:
                if isinstance(renderer, ExternalRenderer):
                    renderer.activate()
        else:
            self._activate_pending_renderers(cls=ExternalRenderer)
            renderers_list = self._default_renderers

        for renderer in renderers_list:
            if isinstance(renderer, ExternalRenderer):
                renderer = copy(renderer)
                for k, v in kwargs.items():
                    if hasattr(renderer, k):
                        setattr(renderer, k, v)

                renderer.render(fig_dict)


# Make renderers a singleton object
# ---------------------------------
renderers = RenderersConfig()
del RenderersConfig


# Show
def show(fig, renderer=None, validate=True, **kwargs):
    """
    Show a figure using either the default renderer(s) or the renderer(s)
    specified by the renderer argument

    Parameters
    ----------
    fig: dict of Figure
        The Figure object or figure dict to display

    renderer: str or None (default None)
        A string containing the names of one or more registered renderers
        (separated by '+' characters) or None.  If None, then the default
        renderers specified in plotly.io.renderers.default are used.

    validate: bool (default True)
        True if the figure should be validated before being shown,
        False otherwise.

    width: int or float
        An integer or float that determines the number of pixels wide the
        plot is. The default is set in plotly.js.

    height: int or float
        An integer or float specifying the height of the plot in pixels.
        The default is set in plotly.js.

    config: dict
        A dict of parameters to configure the figure. The defaults are set
        in plotly.js.

    Returns
    -------
    None
    """
    fig_dict = validate_coerce_fig_to_dict(fig, validate)

    # Mimetype renderers
    bundle = renderers._build_mime_bundle(fig_dict, renderers_string=renderer, **kwargs)
    if bundle:
        if not ipython_display:
            raise ValueError(
                "Mime type rendering requires ipython but it is not installed"
            )

        if not nbformat or Version(nbformat.__version__) < Version("4.2.0"):
            raise ValueError(
                "Mime type rendering requires nbformat>=4.2.0 but it is not installed"
            )

        display_jupyter_version_warnings()

        ipython_display.display(bundle, raw=True)

    # external renderers
    renderers._perform_external_rendering(fig_dict, renderers_string=renderer, **kwargs)


# Register renderers
# ------------------

# Plotly mime type
plotly_renderer = PlotlyRenderer()
renderers["plotly_mimetype"] = plotly_renderer
renderers["jupyterlab"] = plotly_renderer
renderers["nteract"] = plotly_renderer
renderers["vscode"] = plotly_renderer

# HTML-based
config = {}
renderers["notebook"] = NotebookRenderer(config=config)
renderers["notebook_connected"] = NotebookRenderer(config=config, connected=True)
renderers["kaggle"] = KaggleRenderer(config=config)
renderers["azure"] = AzureRenderer(config=config)
renderers["colab"] = ColabRenderer(config=config)
renderers["cocalc"] = CoCalcRenderer()
renderers["databricks"] = DatabricksRenderer()

# JSON
renderers["json"] = JsonRenderer()

# Static Image
renderers["png"] = PngRenderer()
jpeg_renderer = JpegRenderer()
renderers["jpeg"] = jpeg_renderer
renderers["jpg"] = jpeg_renderer
renderers["svg"] = SvgRenderer()
renderers["pdf"] = PdfRenderer()

# External
renderers["browser"] = BrowserRenderer(config=config)
renderers["firefox"] = BrowserRenderer(config=config, using=("firefox"))
renderers["chrome"] = BrowserRenderer(config=config, using=("chrome", "google-chrome"))
renderers["chromium"] = BrowserRenderer(
    config=config, using=("chromium", "chromium-browser")
)
renderers["iframe"] = IFrameRenderer(config=config, include_plotlyjs=True)
renderers["iframe_connected"] = IFrameRenderer(config=config, include_plotlyjs="cdn")
renderers["sphinx_gallery"] = SphinxGalleryHtmlRenderer()
renderers["sphinx_gallery_png"] = SphinxGalleryOrcaRenderer()

# Set default renderer
# --------------------
# Version 4 renderer configuration
default_renderer = None

# Handle the PLOTLY_RENDERER environment variable
env_renderer = os.environ.get("PLOTLY_RENDERER", None)
if env_renderer:
    try:
        renderers._validate_coerce_renderers(env_renderer)
    except ValueError:
        raise ValueError(
            """
Invalid named renderer(s) specified in the 'PLOTLY_RENDERER'
environment variable: {env_renderer}""".format(env_renderer=env_renderer)
        )

    default_renderer = env_renderer
elif ipython and ipython.get_ipython():
    # Try to detect environment so that we can enable a useful
    # default renderer
    if not default_renderer:
        try:
            import google.colab  # noqa: F401

            default_renderer = "colab"
        except ImportError:
            pass

    # Check if we're running in a Kaggle notebook
    if not default_renderer and os.path.exists("/kaggle/input"):
        default_renderer = "kaggle"

    # Check if we're running in an Azure Notebook
    if not default_renderer and "AZURE_NOTEBOOKS_HOST" in os.environ:
        default_renderer = "azure"

    # Check if we're running in VSCode
    if not default_renderer and "VSCODE_PID" in os.environ:
        default_renderer = "vscode"

    # Check if we're running in nteract
    if not default_renderer and "NTERACT_EXE" in os.environ:
        default_renderer = "nteract"

    # Check if we're running in CoCalc
    if not default_renderer and "COCALC_PROJECT_ID" in os.environ:
        default_renderer = "cocalc"

    if not default_renderer and "DATABRICKS_RUNTIME_VERSION" in os.environ:
        default_renderer = "databricks"

    # Check if we're running in spyder and orca is installed
    if not default_renderer and "SPYDER_ARGS" in os.environ:
        try:
            from plotly.io.orca import validate_executable

            validate_executable()
            default_renderer = "svg"
        except ValueError:
            # orca not found
            pass

    # Check if we're running in ipython terminal
    ipython_info = ipython.get_ipython()
    shell = ipython_info.__class__.__name__
    if not default_renderer and (shell == "TerminalInteractiveShell"):
        default_renderer = "browser"

    # Check if we're running in a Jupyter notebook or JupyterLab
    if (
        not default_renderer
        and (shell == "ZMQInteractiveShell")
        and (type(ipython_info).__module__.startswith("ipykernel."))
    ):
        default_renderer = "plotly_mimetype"

    # Fallback to renderer combination that will work automatically
    # in the jupyter notebook, jupyterlab, nteract, vscode, and
    # nbconvert HTML export.
    if not default_renderer:
        default_renderer = "plotly_mimetype+notebook"
else:
    # If ipython isn't available, try to display figures in the default
    # browser
    try:
        import webbrowser

        webbrowser.get()
        default_renderer = "browser"
    except Exception:
        # Many things could have gone wrong
        # There could not be a webbrowser Python module,
        # or the module may be a dumb placeholder
        pass

renderers.render_on_display = True
renderers.default = default_renderer
