from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Fillgradient(_BaseTraceHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "scatter"
    _path_str = "scatter.fillgradient"
    _valid_props = {"colorscale", "start", "stop", "type"}

    # colorscale
    # ----------
    @property
    def colorscale(self):
        """
        Sets the fill gradient colors as a color scale. The color scale
        is interpreted as a gradient applied in the direction specified
        by "orientation", from the lowest to the highest value of the
        scatter plot along that axis, or from the center to the most
        distant point from it, if orientation is "radial".

        The 'colorscale' property is a colorscale and may be
        specified as:
          - A list of colors that will be spaced evenly to create the colorscale.
            Many predefined colorscale lists are included in the sequential, diverging,
            and cyclical modules in the plotly.colors package.
          - A list of 2-element lists where the first element is the
            normalized color level value (starting at 0 and ending at 1),
            and the second item is a valid color string.
            (e.g. [[0, 'green'], [0.5, 'red'], [1.0, 'rgb(0, 0, 255)']])
          - One of the following named colorscales:
                ['aggrnyl', 'agsunset', 'algae', 'amp', 'armyrose', 'balance',
                 'blackbody', 'bluered', 'blues', 'blugrn', 'bluyl', 'brbg',
                 'brwnyl', 'bugn', 'bupu', 'burg', 'burgyl', 'cividis', 'curl',
                 'darkmint', 'deep', 'delta', 'dense', 'earth', 'edge', 'electric',
                 'emrld', 'fall', 'geyser', 'gnbu', 'gray', 'greens', 'greys',
                 'haline', 'hot', 'hsv', 'ice', 'icefire', 'inferno', 'jet',
                 'magenta', 'magma', 'matter', 'mint', 'mrybm', 'mygbm', 'oranges',
                 'orrd', 'oryel', 'oxy', 'peach', 'phase', 'picnic', 'pinkyl',
                 'piyg', 'plasma', 'plotly3', 'portland', 'prgn', 'pubu', 'pubugn',
                 'puor', 'purd', 'purp', 'purples', 'purpor', 'rainbow', 'rdbu',
                 'rdgy', 'rdpu', 'rdylbu', 'rdylgn', 'redor', 'reds', 'solar',
                 'spectral', 'speed', 'sunset', 'sunsetdark', 'teal', 'tealgrn',
                 'tealrose', 'tempo', 'temps', 'thermal', 'tropic', 'turbid',
                 'turbo', 'twilight', 'viridis', 'ylgn', 'ylgnbu', 'ylorbr',
                 'ylorrd'].
            Appending '_r' to a named colorscale reverses it.

        Returns
        -------
        str
        """
        return self["colorscale"]

    @colorscale.setter
    def colorscale(self, val):
        self["colorscale"] = val

    # start
    # -----
    @property
    def start(self):
        """
        Sets the gradient start value. It is given as the absolute
        position on the axis determined by the orientiation. E.g., if
        orientation is "horizontal", the gradient will be horizontal
        and start from the x-position given by start. If omitted, the
        gradient starts at the lowest value of the trace along the
        respective axis. Ignored if orientation is "radial".

        The 'start' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["start"]

    @start.setter
    def start(self, val):
        self["start"] = val

    # stop
    # ----
    @property
    def stop(self):
        """
        Sets the gradient end value. It is given as the absolute
        position on the axis determined by the orientiation. E.g., if
        orientation is "horizontal", the gradient will be horizontal
        and end at the x-position given by end. If omitted, the
        gradient ends at the highest value of the trace along the
        respective axis. Ignored if orientation is "radial".

        The 'stop' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["stop"]

    @stop.setter
    def stop(self, val):
        self["stop"] = val

    # type
    # ----
    @property
    def type(self):
        """
        Sets the type/orientation of the color gradient for the fill.
        Defaults to "none".

        The 'type' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['radial', 'horizontal', 'vertical', 'none']

        Returns
        -------
        Any
        """
        return self["type"]

    @type.setter
    def type(self, val):
        self["type"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        colorscale
            Sets the fill gradient colors as a color scale. The
            color scale is interpreted as a gradient applied in the
            direction specified by "orientation", from the lowest
            to the highest value of the scatter plot along that
            axis, or from the center to the most distant point from
            it, if orientation is "radial".
        start
            Sets the gradient start value. It is given as the
            absolute position on the axis determined by the
            orientiation. E.g., if orientation is "horizontal", the
            gradient will be horizontal and start from the
            x-position given by start. If omitted, the gradient
            starts at the lowest value of the trace along the
            respective axis. Ignored if orientation is "radial".
        stop
            Sets the gradient end value. It is given as the
            absolute position on the axis determined by the
            orientiation. E.g., if orientation is "horizontal", the
            gradient will be horizontal and end at the x-position
            given by end. If omitted, the gradient ends at the
            highest value of the trace along the respective axis.
            Ignored if orientation is "radial".
        type
            Sets the type/orientation of the color gradient for the
            fill. Defaults to "none".
        """

    def __init__(
        self, arg=None, colorscale=None, start=None, stop=None, type=None, **kwargs
    ):
        """
        Construct a new Fillgradient object

        Sets a fill gradient. If not specified, the fillcolor is used
        instead.

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.scatter.Fillgradient`
        colorscale
            Sets the fill gradient colors as a color scale. The
            color scale is interpreted as a gradient applied in the
            direction specified by "orientation", from the lowest
            to the highest value of the scatter plot along that
            axis, or from the center to the most distant point from
            it, if orientation is "radial".
        start
            Sets the gradient start value. It is given as the
            absolute position on the axis determined by the
            orientiation. E.g., if orientation is "horizontal", the
            gradient will be horizontal and start from the
            x-position given by start. If omitted, the gradient
            starts at the lowest value of the trace along the
            respective axis. Ignored if orientation is "radial".
        stop
            Sets the gradient end value. It is given as the
            absolute position on the axis determined by the
            orientiation. E.g., if orientation is "horizontal", the
            gradient will be horizontal and end at the x-position
            given by end. If omitted, the gradient ends at the
            highest value of the trace along the respective axis.
            Ignored if orientation is "radial".
        type
            Sets the type/orientation of the color gradient for the
            fill. Defaults to "none".

        Returns
        -------
        Fillgradient
        """
        super(Fillgradient, self).__init__("fillgradient")

        if "_parent" in kwargs:
            self._parent = kwargs["_parent"]
            return

        # Validate arg
        # ------------
        if arg is None:
            arg = {}
        elif isinstance(arg, self.__class__):
            arg = arg.to_plotly_json()
        elif isinstance(arg, dict):
            arg = _copy.copy(arg)
        else:
            raise ValueError(
                """\
The first argument to the plotly.graph_objs.scatter.Fillgradient
constructor must be a dict or
an instance of :class:`plotly.graph_objs.scatter.Fillgradient`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("colorscale", None)
        _v = colorscale if colorscale is not None else _v
        if _v is not None:
            self["colorscale"] = _v
        _v = arg.pop("start", None)
        _v = start if start is not None else _v
        if _v is not None:
            self["start"] = _v
        _v = arg.pop("stop", None)
        _v = stop if stop is not None else _v
        if _v is not None:
            self["stop"] = _v
        _v = arg.pop("type", None)
        _v = type if type is not None else _v
        if _v is not None:
            self["type"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
