from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Marker(_BaseTraceHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "funnelarea"
    _path_str = "funnelarea.marker"
    _valid_props = {"colors", "colorssrc", "line", "pattern"}

    # colors
    # ------
    @property
    def colors(self):
        """
        Sets the color of each sector. If not specified, the default
        trace color set is used to pick the sector colors.

        The 'colors' property is an array that may be specified as a tuple,
        list, numpy array, or pandas Series

        Returns
        -------
        numpy.ndarray
        """
        return self["colors"]

    @colors.setter
    def colors(self, val):
        self["colors"] = val

    # colorssrc
    # ---------
    @property
    def colorssrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `colors`.

        The 'colorssrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["colorssrc"]

    @colorssrc.setter
    def colorssrc(self, val):
        self["colorssrc"] = val

    # line
    # ----
    @property
    def line(self):
        """
        The 'line' property is an instance of Line
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.funnelarea.marker.Line`
          - A dict of string/value properties that will be passed
            to the Line constructor

            Supported dict properties:

                color
                    Sets the color of the line enclosing each
                    sector. Defaults to the `paper_bgcolor` value.
                colorsrc
                    Sets the source reference on Chart Studio Cloud
                    for `color`.
                width
                    Sets the width (in px) of the line enclosing
                    each sector.
                widthsrc
                    Sets the source reference on Chart Studio Cloud
                    for `width`.

        Returns
        -------
        plotly.graph_objs.funnelarea.marker.Line
        """
        return self["line"]

    @line.setter
    def line(self, val):
        self["line"] = val

    # pattern
    # -------
    @property
    def pattern(self):
        """
        Sets the pattern within the marker.

        The 'pattern' property is an instance of Pattern
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.funnelarea.marker.Pattern`
          - A dict of string/value properties that will be passed
            to the Pattern constructor

            Supported dict properties:

                bgcolor
                    When there is no colorscale sets the color of
                    background pattern fill. Defaults to a
                    `marker.color` background when `fillmode` is
                    "overlay". Otherwise, defaults to a transparent
                    background.
                bgcolorsrc
                    Sets the source reference on Chart Studio Cloud
                    for `bgcolor`.
                fgcolor
                    When there is no colorscale sets the color of
                    foreground pattern fill. Defaults to a
                    `marker.color` background when `fillmode` is
                    "replace". Otherwise, defaults to dark grey or
                    white to increase contrast with the `bgcolor`.
                fgcolorsrc
                    Sets the source reference on Chart Studio Cloud
                    for `fgcolor`.
                fgopacity
                    Sets the opacity of the foreground pattern
                    fill. Defaults to a 0.5 when `fillmode` is
                    "overlay". Otherwise, defaults to 1.
                fillmode
                    Determines whether `marker.color` should be
                    used as a default to `bgcolor` or a `fgcolor`.
                shape
                    Sets the shape of the pattern fill. By default,
                    no pattern is used for filling the area.
                shapesrc
                    Sets the source reference on Chart Studio Cloud
                    for `shape`.
                size
                    Sets the size of unit squares of the pattern
                    fill in pixels, which corresponds to the
                    interval of repetition of the pattern.
                sizesrc
                    Sets the source reference on Chart Studio Cloud
                    for `size`.
                solidity
                    Sets the solidity of the pattern fill. Solidity
                    is roughly the fraction of the area filled by
                    the pattern. Solidity of 0 shows only the
                    background color without pattern and solidty of
                    1 shows only the foreground color without
                    pattern.
                soliditysrc
                    Sets the source reference on Chart Studio Cloud
                    for `solidity`.

        Returns
        -------
        plotly.graph_objs.funnelarea.marker.Pattern
        """
        return self["pattern"]

    @pattern.setter
    def pattern(self, val):
        self["pattern"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        colors
            Sets the color of each sector. If not specified, the
            default trace color set is used to pick the sector
            colors.
        colorssrc
            Sets the source reference on Chart Studio Cloud for
            `colors`.
        line
            :class:`plotly.graph_objects.funnelarea.marker.Line`
            instance or dict with compatible properties
        pattern
            Sets the pattern within the marker.
        """

    def __init__(
        self, arg=None, colors=None, colorssrc=None, line=None, pattern=None, **kwargs
    ):
        """
        Construct a new Marker object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.funnelarea.Marker`
        colors
            Sets the color of each sector. If not specified, the
            default trace color set is used to pick the sector
            colors.
        colorssrc
            Sets the source reference on Chart Studio Cloud for
            `colors`.
        line
            :class:`plotly.graph_objects.funnelarea.marker.Line`
            instance or dict with compatible properties
        pattern
            Sets the pattern within the marker.

        Returns
        -------
        Marker
        """
        super(Marker, self).__init__("marker")

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
The first argument to the plotly.graph_objs.funnelarea.Marker
constructor must be a dict or
an instance of :class:`plotly.graph_objs.funnelarea.Marker`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("colors", None)
        _v = colors if colors is not None else _v
        if _v is not None:
            self["colors"] = _v
        _v = arg.pop("colorssrc", None)
        _v = colorssrc if colorssrc is not None else _v
        if _v is not None:
            self["colorssrc"] = _v
        _v = arg.pop("line", None)
        _v = line if line is not None else _v
        if _v is not None:
            self["line"] = _v
        _v = arg.pop("pattern", None)
        _v = pattern if pattern is not None else _v
        if _v is not None:
            self["pattern"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
