import atexit
import json
import os
import socket
import subprocess
import sys
import threading
import warnings
from copy import copy
from contextlib import contextmanager
from pathlib import Path
from shutil import which

import tenacity

import plotly
from plotly.files import PLOTLY_DIR, ensure_writable_plotly_dir
from plotly.io._utils import validate_coerce_fig_to_dict
from plotly.optional_imports import get_module

psutil = get_module("psutil")

# Valid image format constants
# ----------------------------
valid_formats = ("png", "jpeg", "webp", "svg", "pdf", "eps")
format_conversions = {fmt: fmt for fmt in valid_formats}
format_conversions.update({"jpg": "jpeg"})


# Utility functions
# -----------------
def raise_format_value_error(val):
    raise ValueError(
        """
Invalid value of type {typ} receive as an image format specification.
    Received value: {v}

An image format must be specified as one of the following string values:
    {valid_formats}""".format(
            typ=type(val), v=val, valid_formats=sorted(format_conversions.keys())
        )
    )


def validate_coerce_format(fmt):
    """
    Validate / coerce a user specified image format, and raise an informative
    exception if format is invalid.

    Parameters
    ----------
    fmt
        A value that may or may not be a valid image format string.

    Returns
    -------
    str or None
        A valid image format string as supported by orca. This may not
        be identical to the input image designation. For example,
        the resulting string will always be lower case and  'jpg' is
        converted to 'jpeg'.

        If the input format value is None, then no exception is raised and
        None is returned.

    Raises
    ------
    ValueError
        if the input `fmt` cannot be interpreted as a valid image format.
    """

    # Let None pass through
    if fmt is None:
        return None

    # Check format type
    if not isinstance(fmt, str) or not fmt:
        raise_format_value_error(fmt)

    # Make lower case
    fmt = fmt.lower()

    # Remove leading period, if any.
    # For example '.png' is accepted and converted to 'png'
    if fmt[0] == ".":
        fmt = fmt[1:]

    # Check string value
    if fmt not in format_conversions:
        raise_format_value_error(fmt)

    # Return converted string specification
    return format_conversions[fmt]


def find_open_port():
    """
    Use the socket module to find an open port.

    Returns
    -------
    int
        An open port
    """
    s = socket.socket()
    s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    s.bind(("", 0))
    _, port = s.getsockname()
    s.close()

    return port


# Orca configuration class
# ------------------------
class OrcaConfig(object):
    """
    Singleton object containing the current user defined configuration
    properties for orca.

    These parameters may optionally be saved to the user's ~/.plotly
    directory using the `save` method, in which case they are automatically
    restored in future sessions.
    """

    def __init__(self):

        # Initialize properties dict
        self._props = {}

        # Compute absolute path to the 'plotly/package_data/' directory
        root_dir = os.path.dirname(os.path.abspath(plotly.__file__))
        self.package_dir = os.path.join(root_dir, "package_data")

        # Load pre-existing configuration
        self.reload(warn=False)

        # Compute constants
        plotlyjs = os.path.join(self.package_dir, "plotly.min.js")
        self._constants = {
            "plotlyjs": plotlyjs,
            "config_file": os.path.join(PLOTLY_DIR, ".orca"),
        }

    def restore_defaults(self, reset_server=True):
        """
        Reset all orca configuration properties to their default values
        """
        self._props = {}

        if reset_server:
            # Server must restart before setting is active
            reset_status()

    def update(self, d={}, **kwargs):
        """
        Update one or more properties from a dict or from input keyword
        arguments.

        Parameters
        ----------
        d: dict
            Dictionary from property names to new property values.

        kwargs
            Named argument value pairs where the name is a configuration
            property name and the value is the new property value.

        Returns
        -------
        None

        Examples
        --------
        Update configuration properties using a dictionary

        >>> import plotly.io as pio
        >>> pio.orca.config.update({'timeout': 30, 'default_format': 'svg'})

        Update configuration properties using keyword arguments

        >>> pio.orca.config.update(timeout=30, default_format='svg'})
        """
        # Combine d and kwargs
        if not isinstance(d, dict):
            raise ValueError(
                """
The first argument to update must be a dict, \
but received value of type {typ}l
    Received value: {val}""".format(
                    typ=type(d), val=d
                )
            )

        updates = copy(d)
        updates.update(kwargs)

        # Validate keys
        for k in updates:
            if k not in self._props:
                raise ValueError("Invalid property name: {k}".format(k=k))

        # Apply keys
        for k, v in updates.items():
            setattr(self, k, v)

    def reload(self, warn=True):
        """
        Reload orca settings from ~/.plotly/.orca, if any.

        Note: Settings are loaded automatically when plotly is imported.
        This method is only needed if the setting are changed by some outside
        process (e.g. a text editor) during an interactive session.

        Parameters
        ----------
        warn: bool
            If True, raise informative warnings if settings cannot be restored.
            If False, do not raise warnings if setting cannot be restored.

        Returns
        -------
        None
        """
        if os.path.exists(self.config_file):

            # ### Load file into a string ###
            try:
                with open(self.config_file, "r") as f:
                    orca_str = f.read()
            except:
                if warn:
                    warnings.warn(
                        """\
Unable to read orca configuration file at {path}""".format(
                            path=self.config_file
                        )
                    )
                return

            # ### Parse as JSON ###
            try:
                orca_props = json.loads(orca_str)
            except ValueError:
                if warn:
                    warnings.warn(
                        """\
Orca configuration file at {path} is not valid JSON""".format(
                            path=self.config_file
                        )
                    )
                return

            # ### Update _props ###
            for k, v in orca_props.items():
                self._props[k] = v

        elif warn:
            warnings.warn(
                """\
Orca configuration file at {path} not found""".format(
                    path=self.config_file
                )
            )

    def save(self):
        """
        Attempt to save current settings to disk, so that they are
        automatically restored for future sessions.

        This operation requires write access to the path returned by
        in the `config_file` property.

        Returns
        -------
        None
        """
        if ensure_writable_plotly_dir():
            with open(self.config_file, "w") as f:
                json.dump(self._props, f, indent=4)
        else:
            warnings.warn(
                """\
Failed to write orca configuration file at '{path}'""".format(
                    path=self.config_file
                )
            )

    @property
    def server_url(self):
        """
        The server URL to use for an external orca server, or None if orca
        should be managed locally

        Overrides executable, port, timeout, mathjax, topojson,
        and mapbox_access_token

        Returns
        -------
        str or None
        """
        return self._props.get("server_url", None)

    @server_url.setter
    def server_url(self, val):

        if val is None:
            self._props.pop("server_url", None)
            return
        if not isinstance(val, str):
            raise ValueError(
                """
The server_url property must be a string, but received value of type {typ}.
    Received value: {val}""".format(
                    typ=type(val), val=val
                )
            )

        if not val.startswith("http://") and not val.startswith("https://"):
            val = "http://" + val

        shutdown_server()
        self.executable = None
        self.port = None
        self.timeout = None
        self.mathjax = None
        self.topojson = None
        self.mapbox_access_token = None
        self._props["server_url"] = val

    @property
    def port(self):
        """
        The specific port to use to communicate with the orca server, or
        None if the port is to be chosen automatically.

        If an orca server is active, the port in use is stored in the
        plotly.io.orca.status.port property.

        Returns
        -------
        int or None
        """
        return self._props.get("port", None)

    @port.setter
    def port(self, val):

        if val is None:
            self._props.pop("port", None)
            return
        if not isinstance(val, int):
            raise ValueError(
                """
The port property must be an integer, but received value of type {typ}.
    Received value: {val}""".format(
                    typ=type(val), val=val
                )
            )

        self._props["port"] = val

    @property
    def executable(self):
        """
        The name or full path of the orca executable.

         - If a name (e.g. 'orca'), then it should be the name of an orca
           executable on the PATH. The directories on the PATH can be
           displayed by running the following command:

           >>> import os
           >>> print(os.environ.get('PATH').replace(os.pathsep, os.linesep))

         - If a full path (e.g. '/path/to/orca'), then
           it should be the full path to an orca executable. In this case
           the executable does not need to reside on the PATH.

        If an orca server has been validated, then the full path to the
        validated orca executable is stored in the
        plotly.io.orca.status.executable property.

        Returns
        -------
        str
        """
        executable_list = self._props.get("executable_list", ["orca"])
        if executable_list is None:
            return None
        else:
            return " ".join(executable_list)

    @executable.setter
    def executable(self, val):

        if val is None:
            self._props.pop("executable", None)
        else:
            if not isinstance(val, str):
                raise ValueError(
                    """
The executable property must be a string, but received value of type {typ}.
    Received value: {val}""".format(
                        typ=type(val), val=val
                    )
                )
            if isinstance(val, str):
                val = [val]
            self._props["executable_list"] = val

        # Server and validation must restart before setting is active
        reset_status()

    @property
    def timeout(self):
        """
        The number of seconds of inactivity required before the orca server
        is shut down.

        For example, if timeout is set to 20, then the orca
        server will shutdown once is has not been used for at least
        20 seconds. If timeout is set to None, then the server will not be
        automatically shut down due to inactivity.

        Regardless of the value of timeout, a running orca server may be
        manually shut down like this:

        >>> import plotly.io as pio
        >>> pio.orca.shutdown_server()

        Returns
        -------
        int or float or None
        """
        return self._props.get("timeout", None)

    @timeout.setter
    def timeout(self, val):

        if val is None:
            self._props.pop("timeout", None)
        else:
            if not isinstance(val, (int, float)):
                raise ValueError(
                    """
The timeout property must be a number, but received value of type {typ}.
    Received value: {val}""".format(
                        typ=type(val), val=val
                    )
                )
            self._props["timeout"] = val

        # Server must restart before setting is active
        shutdown_server()

    @property
    def default_width(self):
        """
        The default width to use on image export. This value is only
        applied if no width value is supplied to the plotly.io
        to_image or write_image functions.

        Returns
        -------
        int or None
        """
        return self._props.get("default_width", None)

    @default_width.setter
    def default_width(self, val):

        if val is None:
            self._props.pop("default_width", None)
            return
        if not isinstance(val, int):
            raise ValueError(
                """
The default_width property must be an int, but received value of type {typ}.
    Received value: {val}""".format(
                    typ=type(val), val=val
                )
            )
        self._props["default_width"] = val

    @property
    def default_height(self):
        """
        The default height to use on image export. This value is only
        applied if no height value is supplied to the plotly.io
        to_image or write_image functions.

        Returns
        -------
        int or None
        """
        return self._props.get("default_height", None)

    @default_height.setter
    def default_height(self, val):

        if val is None:
            self._props.pop("default_height", None)
            return
        if not isinstance(val, int):
            raise ValueError(
                """
The default_height property must be an int, but received value of type {typ}.
    Received value: {val}""".format(
                    typ=type(val), val=val
                )
            )
        self._props["default_height"] = val

    @property
    def default_format(self):
        """
        The default image format to use on image export.

        Valid image formats strings are:
          - 'png'
          - 'jpg' or 'jpeg'
          - 'webp'
          - 'svg'
          - 'pdf'
          - 'eps' (Requires the poppler library to be installed)

        This value is only applied if no format value is supplied to the
        plotly.io to_image or write_image functions.

        Returns
        -------
        str or None
        """
        return self._props.get("default_format", "png")

    @default_format.setter
    def default_format(self, val):
        if val is None:
            self._props.pop("default_format", None)
            return

        val = validate_coerce_format(val)
        self._props["default_format"] = val

    @property
    def default_scale(self):
        """
        The default image scaling factor to use on image export.
        This value is only applied if no scale value is supplied to the
        plotly.io to_image or write_image functions.

        Returns
        -------
        int or None
        """
        return self._props.get("default_scale", 1)

    @default_scale.setter
    def default_scale(self, val):

        if val is None:
            self._props.pop("default_scale", None)
            return
        if not isinstance(val, (int, float)):
            raise ValueError(
                """
The default_scale property must be a number, but received value of type {typ}.
    Received value: {val}""".format(
                    typ=type(val), val=val
                )
            )
        self._props["default_scale"] = val

    @property
    def topojson(self):
        """
        Path to the topojson files needed to render choropleth traces.

        If None, topojson files from the plot.ly CDN are used.

        Returns
        -------
        str
        """
        return self._props.get("topojson", None)

    @topojson.setter
    def topojson(self, val):

        if val is None:
            self._props.pop("topojson", None)
        else:
            if not isinstance(val, str):
                raise ValueError(
                    """
The topojson property must be a string, but received value of type {typ}.
    Received value: {val}""".format(
                        typ=type(val), val=val
                    )
                )
            self._props["topojson"] = val

        # Server must restart before setting is active
        shutdown_server()

    @property
    def mathjax(self):
        """
        Path to the MathJax bundle needed to render LaTeX characters

        Returns
        -------
        str
        """
        return self._props.get(
            "mathjax",
            ("https://cdnjs.cloudflare.com" "/ajax/libs/mathjax/2.7.5/MathJax.js"),
        )

    @mathjax.setter
    def mathjax(self, val):

        if val is None:
            self._props.pop("mathjax", None)
        else:
            if not isinstance(val, str):
                raise ValueError(
                    """
The mathjax property must be a string, but received value of type {typ}.
    Received value: {val}""".format(
                        typ=type(val), val=val
                    )
                )
            self._props["mathjax"] = val

        # Server must restart before setting is active
        shutdown_server()

    @property
    def mapbox_access_token(self):
        """
        Mapbox access token required to render mapbox traces.

        Returns
        -------
        str
        """
        return self._props.get("mapbox_access_token", None)

    @mapbox_access_token.setter
    def mapbox_access_token(self, val):

        if val is None:
            self._props.pop("mapbox_access_token", None)
        else:
            if not isinstance(val, str):
                raise ValueError(
                    """
The mapbox_access_token property must be a string, \
but received value of type {typ}.
    Received value: {val}""".format(
                        typ=type(val), val=val
                    )
                )
            self._props["mapbox_access_token"] = val

        # Server must restart before setting is active
        shutdown_server()

    @property
    def use_xvfb(self):
        dflt = "auto"
        return self._props.get("use_xvfb", dflt)

    @use_xvfb.setter
    def use_xvfb(self, val):
        valid_vals = [True, False, "auto"]
        if val is None:
            self._props.pop("use_xvfb", None)
        else:
            if val not in valid_vals:
                raise ValueError(
                    """
The use_xvfb property must be one of {valid_vals}
    Received value of type {typ}: {val}""".format(
                        valid_vals=valid_vals, typ=type(val), val=repr(val)
                    )
                )

            self._props["use_xvfb"] = val

        # Server and validation must restart before setting is active
        reset_status()

    @property
    def plotlyjs(self):
        """
        The plotly.js bundle being used for image rendering.

        Returns
        -------
        str
        """
        return self._constants.get("plotlyjs", None)

    @property
    def config_file(self):
        """
        Path to orca configuration file

        Using the `plotly.io.config.save()` method will save the current
        configuration settings to this file. Settings in this file are
        restored at the beginning of each sessions.

        Returns
        -------
        str
        """
        return os.path.join(PLOTLY_DIR, ".orca")

    def __repr__(self):
        """
        Display a nice representation of the current orca configuration.
        """
        return """\
orca configuration
------------------
    server_url: {server_url}
    executable: {executable}
    port: {port}
    timeout: {timeout}
    default_width: {default_width}
    default_height: {default_height}
    default_scale: {default_scale}
    default_format: {default_format}
    mathjax: {mathjax}
    topojson: {topojson}
    mapbox_access_token: {mapbox_access_token}
    use_xvfb: {use_xvfb}

constants
---------
    plotlyjs: {plotlyjs}
    config_file: {config_file}

""".format(
            server_url=self.server_url,
            port=self.port,
            executable=self.executable,
            timeout=self.timeout,
            default_width=self.default_width,
            default_height=self.default_height,
            default_scale=self.default_scale,
            default_format=self.default_format,
            mathjax=self.mathjax,
            topojson=self.topojson,
            mapbox_access_token=self.mapbox_access_token,
            plotlyjs=self.plotlyjs,
            config_file=self.config_file,
            use_xvfb=self.use_xvfb,
        )


# Make config a singleton object
# ------------------------------
config = OrcaConfig()
del OrcaConfig


# Orca status class
# ------------------------
class OrcaStatus(object):
    """
    Class to store information about the current status of the orca server.
    """

    _props = {
        "state": "unvalidated",  # or 'validated' or 'running'
        "executable_list": None,
        "version": None,
        "pid": None,
        "port": None,
        "command": None,
    }

    @property
    def state(self):
        """
        A string representing the state of the orca server process

        One of:
          - unvalidated: The orca executable has not yet been searched for or
            tested to make sure its valid.
          - validated: The orca executable has been located and tested for
            validity, but it is not running.
          - running: The orca server process is currently running.
        """
        return self._props["state"]

    @property
    def executable(self):
        """
        If the `state` property is 'validated' or 'running', this property
        contains the full path to the orca executable.

        This path can be specified explicitly by setting the `executable`
        property of the `plotly.io.orca.config` object.

        This property will be None if the `state` is 'unvalidated'.
        """
        executable_list = self._props["executable_list"]
        if executable_list is None:
            return None
        else:
            return " ".join(executable_list)

    @property
    def version(self):
        """
        If the `state` property is 'validated' or 'running', this property
        contains the version of the validated orca executable.

        This property will be None if the `state` is 'unvalidated'.
        """
        return self._props["version"]

    @property
    def pid(self):
        """
        The process id of the orca server process, if any. This property
        will be None if the `state` is not 'running'.
        """
        return self._props["pid"]

    @property
    def port(self):
        """
        The port number that the orca server process is listening to, if any.
        This property will be None if the `state` is not 'running'.

        This port can be specified explicitly by setting the `port`
        property of the `plotly.io.orca.config` object.
        """
        return self._props["port"]

    @property
    def command(self):
        """
        The command arguments used to launch the running orca server, if any.
        This property will be None if the `state` is not 'running'.
        """
        return self._props["command"]

    def __repr__(self):
        """
        Display a nice representation of the current orca server status.
        """
        return """\
orca status
-----------
    state: {state}
    executable: {executable}
    version: {version}
    port: {port}
    pid: {pid}
    command: {command}

""".format(
            executable=self.executable,
            version=self.version,
            port=self.port,
            pid=self.pid,
            state=self.state,
            command=self.command,
        )


# Make status a singleton object
# ------------------------------
status = OrcaStatus()
del OrcaStatus


@contextmanager
def orca_env():
    """
    Context manager to clear and restore environment variables that are
    problematic for orca to function properly

    NODE_OPTIONS: When this variable is set, orca <v1.2 will have a
    segmentation fault due to an electron bug.
    See: https://github.com/electron/electron/issues/12695

    ELECTRON_RUN_AS_NODE: When this environment variable is set the call
    to orca is transformed into a call to nodejs.
    See https://github.com/plotly/orca/issues/149#issuecomment-443506732
    """
    clear_env_vars = ["NODE_OPTIONS", "ELECTRON_RUN_AS_NODE", "LD_PRELOAD"]
    orig_env_vars = {}

    try:
        # Clear and save
        orig_env_vars.update(
            {var: os.environ.pop(var) for var in clear_env_vars if var in os.environ}
        )
        yield
    finally:
        # Restore
        for var, val in orig_env_vars.items():
            os.environ[var] = val


# Public orca server interaction functions
# ----------------------------------------
def validate_executable():
    """
    Attempt to find and validate the orca executable specified by the
    `plotly.io.orca.config.executable` property.

    If the `plotly.io.orca.status.state` property is 'validated' or 'running'
    then this function does nothing.

    How it works:
      - First, it searches the system PATH for an executable that matches the
      name or path specified in the `plotly.io.orca.config.executable`
      property.
      - Then it runs the executable with the `--help` flag to make sure
      it's the plotly orca executable
      - Then it runs the executable with the `--version` flag to check the
      orca version.

    If all of these steps are successful then the `status.state` property
    is set to 'validated' and the `status.executable` and `status.version`
    properties are populated

    Returns
    -------
    None
    """
    # Check state
    # -----------
    if status.state != "unvalidated":
        # Nothing more to do
        return

    # Initialize error messages
    # -------------------------
    install_location_instructions = """\
If you haven't installed orca yet, you can do so using conda as follows:

    $ conda install -c plotly plotly-orca

Alternatively, see other installation methods in the orca project README at
https://github.com/plotly/orca

After installation is complete, no further configuration should be needed.

If you have installed orca, then for some reason plotly.py was unable to
locate it. In this case, set the `plotly.io.orca.config.executable`
property to the full path of your orca executable. For example:

    >>> plotly.io.orca.config.executable = '/path/to/orca'

After updating this executable property, try the export operation again.
If it is successful then you may want to save this configuration so that it
will be applied automatically in future sessions. You can do this as follows:

    >>> plotly.io.orca.config.save()

If you're still having trouble, feel free to ask for help on the forums at
https://community.plot.ly/c/api/python
"""

    # Try to find an executable
    # -------------------------
    # Search for executable name or path in config.executable
    executable = which(config.executable)
    path = os.environ.get("PATH", os.defpath)
    formatted_path = path.replace(os.pathsep, "\n    ")

    if executable is None:
        raise ValueError(
            """
The orca executable is required to export figures as static images,
but it could not be found on the system path.

Searched for executable '{executable}' on the following path:
    {formatted_path}

{instructions}""".format(
                executable=config.executable,
                formatted_path=formatted_path,
                instructions=install_location_instructions,
            )
        )

    # Check if we should run with Xvfb
    # --------------------------------
    xvfb_args = [
        "--auto-servernum",
        "--server-args",
        "-screen 0 640x480x24 +extension RANDR +extension GLX",
        executable,
    ]

    if config.use_xvfb == True:
        # Use xvfb
        xvfb_run_executable = which("xvfb-run")
        if not xvfb_run_executable:
            raise ValueError(
                """
The plotly.io.orca.config.use_xvfb property is set to True, but the
xvfb-run executable could not be found on the system path.

Searched for the executable 'xvfb-run' on the following path:
    {formatted_path}""".format(
                    formatted_path=formatted_path
                )
            )

        executable_list = [xvfb_run_executable] + xvfb_args
    elif (
        config.use_xvfb == "auto"
        and sys.platform.startswith("linux")
        and not os.environ.get("DISPLAY")
        and which("xvfb-run")
    ):
        # use_xvfb is 'auto', we're on linux without a display server,
        # and xvfb-run is available. Use it.
        xvfb_run_executable = which("xvfb-run")
        executable_list = [xvfb_run_executable] + xvfb_args
    else:
        # Do not use xvfb
        executable_list = [executable]

    # Run executable with --help and see if it's our orca
    # ---------------------------------------------------
    invalid_executable_msg = """
The orca executable is required in order to export figures as static images,
but the executable that was found at '{executable}'
does not seem to be a valid plotly orca executable. Please refer to the end of
this message for details on what went wrong.

{instructions}""".format(
        executable=executable, instructions=install_location_instructions
    )

    # ### Run with Popen so we get access to stdout and stderr
    with orca_env():
        p = subprocess.Popen(
            executable_list + ["--help"], stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )

        help_result, help_error = p.communicate()

    if p.returncode != 0:
        err_msg = (
            invalid_executable_msg
            + """
Here is the error that was returned by the command
    $ {executable} --help

[Return code: {returncode}]
{err_msg}
""".format(
                executable=" ".join(executable_list),
                err_msg=help_error.decode("utf-8"),
                returncode=p.returncode,
            )
        )

        # Check for Linux without X installed.
        if sys.platform.startswith("linux") and not os.environ.get("DISPLAY"):

            err_msg += """\
Note: When used on Linux, orca requires an X11 display server, but none was
detected. Please install Xvfb and configure plotly.py to run orca using Xvfb
as follows:

    >>> import plotly.io as pio
    >>> pio.orca.config.use_xvfb = True

You can save this configuration for use in future sessions as follows:

    >>> pio.orca.config.save()

See https://www.x.org/releases/X11R7.6/doc/man/man1/Xvfb.1.xhtml
for more info on Xvfb
"""
        raise ValueError(err_msg)

    if not help_result:
        raise ValueError(
            invalid_executable_msg
            + """
The error encountered is that no output was returned by the command
    $ {executable} --help
""".format(
                executable=" ".join(executable_list)
            )
        )

    if "Plotly's image-exporting utilities" not in help_result.decode("utf-8"):
        raise ValueError(
            invalid_executable_msg
            + """
The error encountered is that unexpected output was returned by the command
    $ {executable} --help

{help_result}
""".format(
                executable=" ".join(executable_list), help_result=help_result
            )
        )

    # Get orca version
    # ----------------
    # ### Run with Popen so we get access to stdout and stderr
    with orca_env():
        p = subprocess.Popen(
            executable_list + ["--version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        version_result, version_error = p.communicate()

    if p.returncode != 0:
        raise ValueError(
            invalid_executable_msg
            + """
An error occurred while trying to get the version of the orca executable.
Here is the command that plotly.py ran to request the version
    $ {executable} --version

This command returned the following error:

[Return code: {returncode}]
{err_msg}
        """.format(
                executable=" ".join(executable_list),
                err_msg=version_error.decode("utf-8"),
                returncode=p.returncode,
            )
        )

    if not version_result:
        raise ValueError(
            invalid_executable_msg
            + """
The error encountered is that no version was reported by the orca executable.
Here is the command that plotly.py ran to request the version:

    $ {executable} --version
""".format(
                executable=" ".join(executable_list)
            )
        )
    else:
        version_result = version_result.decode()

    status._props["executable_list"] = executable_list
    status._props["version"] = version_result.strip()
    status._props["state"] = "validated"


def reset_status():
    """
    Shutdown the running orca server, if any, and reset the orca status
    to unvalidated.

    This command is only needed if the desired orca executable is changed
    during an interactive session.

    Returns
    -------
    None
    """
    shutdown_server()
    status._props["executable_list"] = None
    status._props["version"] = None
    status._props["state"] = "unvalidated"


# Initialze process control variables
# -----------------------------------
orca_lock = threading.Lock()
orca_state = {"proc": None, "shutdown_timer": None}


# Shutdown
# --------
# The @atexit.register annotation ensures that the shutdown function is
# is run when the Python process is terminated
@atexit.register
def cleanup():
    shutdown_server()


def shutdown_server():
    """
    Shutdown the running orca server process, if any

    Returns
    -------
    None
    """
    # Use double-check locking to make sure the properties of orca_state
    # are updated consistently across threads.
    if orca_state["proc"] is not None:
        with orca_lock:
            if orca_state["proc"] is not None:

                # We use psutil to kill all child processes of the main orca
                # process. This prevents any zombie processes from being
                # left over, and it saves us from needing to write
                # OS-specific process management code here.

                parent = psutil.Process(orca_state["proc"].pid)
                for child in parent.children(recursive=True):
                    try:
                        child.terminate()
                    except:
                        # We tried, move on
                        pass

                try:
                    # Kill parent process
                    orca_state["proc"].terminate()

                    # Wait for the process to shutdown
                    child_status = orca_state["proc"].wait()
                except:
                    # We tried, move on
                    pass

                # Update our internal process management state
                orca_state["proc"] = None

                if orca_state["shutdown_timer"] is not None:
                    orca_state["shutdown_timer"].cancel()
                    orca_state["shutdown_timer"] = None

                orca_state["port"] = None

                # Update orca.status so the user has an accurate view
                # of the state of the orca server
                status._props["state"] = "validated"
                status._props["pid"] = None
                status._props["port"] = None
                status._props["command"] = None


# Launch or get server
def ensure_server():
    """
    Start an orca server if none is running. If a server is already running,
    then reset the timeout countdown

    Returns
    -------
    None
    """

    # Validate psutil
    if psutil is None:
        raise ValueError(
            """\
Image generation requires the psutil package.

Install using pip:
    $ pip install psutil

Install using conda:
    $ conda install psutil
"""
        )

    # Validate requests
    if not get_module("requests"):
        raise ValueError(
            """\
Image generation requires the requests package.

Install using pip:
    $ pip install requests

Install using conda:
    $ conda install requests
"""
        )

    if not config.server_url:
        # Validate orca executable only if server_url is not provided
        if status.state == "unvalidated":
            validate_executable()
        # Acquire lock to make sure that we keep the properties of orca_state
        # consistent across threads
        with orca_lock:
            # Cancel the current shutdown timer, if any
            if orca_state["shutdown_timer"] is not None:
                orca_state["shutdown_timer"].cancel()

            # Start a new server process if none is active
            if orca_state["proc"] is None:

                # Determine server port
                if config.port is None:
                    orca_state["port"] = find_open_port()
                else:
                    orca_state["port"] = config.port

                # Build orca command list
                cmd_list = status._props["executable_list"] + [
                    "serve",
                    "-p",
                    str(orca_state["port"]),
                    "--plotly",
                    config.plotlyjs,
                    "--graph-only",
                ]

                if config.topojson:
                    cmd_list.extend(["--topojson", config.topojson])

                if config.mathjax:
                    cmd_list.extend(["--mathjax", config.mathjax])

                if config.mapbox_access_token:
                    cmd_list.extend(
                        ["--mapbox-access-token", config.mapbox_access_token]
                    )

                # Create subprocess that launches the orca server on the
                # specified port.
                DEVNULL = open(os.devnull, "wb")
                with orca_env():
                    stderr = DEVNULL if "CI" in os.environ else None  # fix for CI
                    orca_state["proc"] = subprocess.Popen(
                        cmd_list, stdout=DEVNULL, stderr=stderr
                    )

                # Update orca.status so the user has an accurate view
                # of the state of the orca server
                status._props["state"] = "running"
                status._props["pid"] = orca_state["proc"].pid
                status._props["port"] = orca_state["port"]
                status._props["command"] = cmd_list

            # Create new shutdown timer if a timeout was specified
            if config.timeout is not None:
                t = threading.Timer(config.timeout, shutdown_server)
                # Make it a daemon thread so that exit won't wait for timer to
                # complete
                t.daemon = True
                t.start()
                orca_state["shutdown_timer"] = t


@tenacity.retry(
    wait=tenacity.wait_random(min=5, max=10),
    stop=tenacity.stop_after_delay(60000),
)
def request_image_with_retrying(**kwargs):
    """
    Helper method to perform an image request to a running orca server process
    with retrying logic.
    """
    from requests import post
    from plotly.io.json import to_json_plotly

    if config.server_url:
        server_url = config.server_url
    else:
        server_url = "http://{hostname}:{port}".format(
            hostname="localhost", port=orca_state["port"]
        )

    request_params = {k: v for k, v, in kwargs.items() if v is not None}
    json_str = to_json_plotly(request_params)
    response = post(server_url + "/", data=json_str)

    if response.status_code == 522:
        # On "522: client socket timeout", return server and keep trying
        shutdown_server()
        ensure_server()
        raise OSError("522: client socket timeout")

    return response


def to_image(fig, format=None, width=None, height=None, scale=None, validate=True):
    """
    Convert a figure to a static image bytes string

    Parameters
    ----------
    fig:
        Figure object or dict representing a figure

    format: str or None
        The desired image format. One of
          - 'png'
          - 'jpg' or 'jpeg'
          - 'webp'
          - 'svg'
          - 'pdf'
          - 'eps' (Requires the poppler library to be installed)

        If not specified, will default to `plotly.io.config.default_format`

    width: int or None
        The width of the exported image in layout pixels. If the `scale`
        property is 1.0, this will also be the width of the exported image
        in physical pixels.

        If not specified, will default to `plotly.io.config.default_width`

    height: int or None
        The height of the exported image in layout pixels. If the `scale`
        property is 1.0, this will also be the height of the exported image
        in physical pixels.

        If not specified, will default to `plotly.io.config.default_height`

    scale: int or float or None
        The scale factor to use when exporting the figure. A scale factor
        larger than 1.0 will increase the image resolution with respect
        to the figure's layout pixel dimensions. Whereas as scale factor of
        less than 1.0 will decrease the image resolution.

        If not specified, will default to `plotly.io.config.default_scale`

    validate: bool
        True if the figure should be validated before being converted to
        an image, False otherwise.

    Returns
    -------
    bytes
        The image data
    """
    # Make sure orca sever is running
    # -------------------------------
    ensure_server()

    # Handle defaults
    # ---------------
    # Apply configuration defaults to unspecified arguments
    if format is None:
        format = config.default_format

    format = validate_coerce_format(format)

    if scale is None:
        scale = config.default_scale

    if width is None:
        width = config.default_width

    if height is None:
        height = config.default_height

    # Validate figure
    # ---------------
    fig_dict = validate_coerce_fig_to_dict(fig, validate)

    # Request image from server
    # -------------------------
    try:
        response = request_image_with_retrying(
            figure=fig_dict, format=format, scale=scale, width=width, height=height
        )
    except OSError as err:
        # Get current status string
        status_str = repr(status)

        if config.server_url:
            raise ValueError(
                """
Plotly.py was unable to communicate with the orca server at {server_url}

Please check that the server is running and accessible.
""".format(
                    server_url=config.server_url
                )
            )

        else:
            # Check if the orca server process exists
            pid_exists = psutil.pid_exists(status.pid)

            # Raise error message based on whether the server process existed
            if pid_exists:
                raise ValueError(
                    """
For some reason plotly.py was unable to communicate with the
local orca server process, even though the server process seems to be running.

Please review the process and connection information below:

{info}
""".format(
                        info=status_str
                    )
                )
            else:
                # Reset the status so that if the user tries again, we'll try to
                # start the server again
                reset_status()
                raise ValueError(
                    """
For some reason the orca server process is no longer running.

Please review the process and connection information below:

{info}
plotly.py will attempt to start the local server process again the next time
an image export operation is performed.
""".format(
                        info=status_str
                    )
                )

    # Check response
    # --------------
    if response.status_code == 200:
        # All good
        return response.content
    else:
        # ### Something went wrong ###
        err_message = """
The image request was rejected by the orca conversion utility
with the following error:
   {status}: {msg}
""".format(
            status=response.status_code, msg=response.content.decode("utf-8")
        )

        # ### Try to be helpful ###
        # Status codes from /src/component/plotly-graph/constants.js in the
        # orca code base.
        # statusMsg: {
        #     400: 'invalid or malformed request syntax',
        #     522: client socket timeout
        #     525: 'plotly.js error',
        #     526: 'plotly.js version 1.11.0 or up required',
        #     530: 'image conversion error'
        # }
        if response.status_code == 400 and isinstance(fig, dict) and not validate:
            err_message += """
Try setting the `validate` argument to True to check for errors in the
figure specification"""
        elif response.status_code == 525:
            any_mapbox = any(
                [
                    trace.get("type", None) == "scattermapbox"
                    for trace in fig_dict.get("data", [])
                ]
            )
            if any_mapbox and config.mapbox_access_token is None:
                err_message += """
Exporting scattermapbox traces requires a mapbox access token.
Create a token in your mapbox account and then set it using:

>>> plotly.io.orca.config.mapbox_access_token = 'pk.abc...'

If you would like this token to be applied automatically in
future sessions, then save your orca configuration as follows:

>>> plotly.io.orca.config.save()
"""
        elif response.status_code == 530 and format == "eps":
            err_message += """
Exporting to EPS format requires the poppler library.  You can install
poppler on MacOS or Linux with:

    $ conda install poppler

Or, you can install it on MacOS using homebrew with:

    $ brew install poppler

Or, you can install it on Linux using your distribution's package manager to
install the 'poppler-utils' package.

Unfortunately, we don't yet know of an easy way to install poppler on Windows.
"""
        raise ValueError(err_message)


def write_image(
    fig, file, format=None, scale=None, width=None, height=None, validate=True
):
    """
    Convert a figure to a static image and write it to a file or writeable
    object

    Parameters
    ----------
    fig:
        Figure object or dict representing a figure

    file: str or writeable
        A string representing a local file path or a writeable object
        (e.g. a pathlib.Path object or an open file descriptor)

    format: str or None
        The desired image format. One of
          - 'png'
          - 'jpg' or 'jpeg'
          - 'webp'
          - 'svg'
          - 'pdf'
          - 'eps' (Requires the poppler library to be installed)

        If not specified and `file` is a string then this will default to the
        file extension. If not specified and `file` is not a string then this
        will default to `plotly.io.config.default_format`

    width: int or None
        The width of the exported image in layout pixels. If the `scale`
        property is 1.0, this will also be the width of the exported image
        in physical pixels.

        If not specified, will default to `plotly.io.config.default_width`

    height: int or None
        The height of the exported image in layout pixels. If the `scale`
        property is 1.0, this will also be the height of the exported image
        in physical pixels.

        If not specified, will default to `plotly.io.config.default_height`

    scale: int or float or None
        The scale factor to use when exporting the figure. A scale factor
        larger than 1.0 will increase the image resolution with respect
        to the figure's layout pixel dimensions. Whereas as scale factor of
        less than 1.0 will decrease the image resolution.

        If not specified, will default to `plotly.io.config.default_scale`

    validate: bool
        True if the figure should be validated before being converted to
        an image, False otherwise.

    Returns
    -------
    None
    """

    # Try to cast `file` as a pathlib object `path`.
    # ----------------------------------------------
    if isinstance(file, str):
        # Use the standard Path constructor to make a pathlib object.
        path = Path(file)
    elif isinstance(file, Path):
        # `file` is already a Path object.
        path = file
    else:
        # We could not make a Path object out of file. Either `file` is an open file
        # descriptor with a `write()` method or it's an invalid object.
        path = None

    # Infer format if not specified
    # -----------------------------
    if path is not None and format is None:
        ext = path.suffix
        if ext:
            format = ext.lstrip(".")
        else:
            raise ValueError(
                """
Cannot infer image type from output path '{file}'.
Please add a file extension or specify the type using the format parameter.
For example:

    >>> import plotly.io as pio
    >>> pio.write_image(fig, file_path, format='png')
""".format(
                    file=file
                )
            )

    # Request image
    # -------------
    # Do this first so we don't create a file if image conversion fails
    img_data = to_image(
        fig, format=format, scale=scale, width=width, height=height, validate=validate
    )

    # Open file
    # ---------
    if path is None:
        # We previously failed to make sense of `file` as a pathlib object.
        # Attempt to write to `file` as an open file descriptor.
        try:
            file.write(img_data)
            return
        except AttributeError:
            pass
        raise ValueError(
            """
The 'file' argument '{file}' is not a string, pathlib.Path object, or file descriptor.
""".format(
                file=file
            )
        )
    else:
        # We previously succeeded in interpreting `file` as a pathlib object.
        # Now we can use `write_bytes()`.
        path.write_bytes(img_data)
