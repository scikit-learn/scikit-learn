from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Contours(_BaseTraceHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "surface"
    _path_str = "surface.contours"
    _valid_props = {"x", "y", "z"}

    # x
    # -
    @property
    def x(self):
        """
        The 'x' property is an instance of X
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.surface.contours.X`
          - A dict of string/value properties that will be passed
            to the X constructor

            Supported dict properties:

                color
                    Sets the color of the contour lines.
                end
                    Sets the end contour level value. Must be more
                    than `contours.start`
                highlight
                    Determines whether or not contour lines about
                    the x dimension are highlighted on hover.
                highlightcolor
                    Sets the color of the highlighted contour
                    lines.
                highlightwidth
                    Sets the width of the highlighted contour
                    lines.
                project
                    :class:`plotly.graph_objects.surface.contours.x
                    .Project` instance or dict with compatible
                    properties
                show
                    Determines whether or not contour lines about
                    the x dimension are drawn.
                size
                    Sets the step between each contour level. Must
                    be positive.
                start
                    Sets the starting contour level value. Must be
                    less than `contours.end`
                usecolormap
                    An alternate to "color". Determines whether or
                    not the contour lines are colored using the
                    trace "colorscale".
                width
                    Sets the width of the contour lines.

        Returns
        -------
        plotly.graph_objs.surface.contours.X
        """
        return self["x"]

    @x.setter
    def x(self, val):
        self["x"] = val

    # y
    # -
    @property
    def y(self):
        """
        The 'y' property is an instance of Y
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.surface.contours.Y`
          - A dict of string/value properties that will be passed
            to the Y constructor

            Supported dict properties:

                color
                    Sets the color of the contour lines.
                end
                    Sets the end contour level value. Must be more
                    than `contours.start`
                highlight
                    Determines whether or not contour lines about
                    the y dimension are highlighted on hover.
                highlightcolor
                    Sets the color of the highlighted contour
                    lines.
                highlightwidth
                    Sets the width of the highlighted contour
                    lines.
                project
                    :class:`plotly.graph_objects.surface.contours.y
                    .Project` instance or dict with compatible
                    properties
                show
                    Determines whether or not contour lines about
                    the y dimension are drawn.
                size
                    Sets the step between each contour level. Must
                    be positive.
                start
                    Sets the starting contour level value. Must be
                    less than `contours.end`
                usecolormap
                    An alternate to "color". Determines whether or
                    not the contour lines are colored using the
                    trace "colorscale".
                width
                    Sets the width of the contour lines.

        Returns
        -------
        plotly.graph_objs.surface.contours.Y
        """
        return self["y"]

    @y.setter
    def y(self, val):
        self["y"] = val

    # z
    # -
    @property
    def z(self):
        """
        The 'z' property is an instance of Z
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.surface.contours.Z`
          - A dict of string/value properties that will be passed
            to the Z constructor

            Supported dict properties:

                color
                    Sets the color of the contour lines.
                end
                    Sets the end contour level value. Must be more
                    than `contours.start`
                highlight
                    Determines whether or not contour lines about
                    the z dimension are highlighted on hover.
                highlightcolor
                    Sets the color of the highlighted contour
                    lines.
                highlightwidth
                    Sets the width of the highlighted contour
                    lines.
                project
                    :class:`plotly.graph_objects.surface.contours.z
                    .Project` instance or dict with compatible
                    properties
                show
                    Determines whether or not contour lines about
                    the z dimension are drawn.
                size
                    Sets the step between each contour level. Must
                    be positive.
                start
                    Sets the starting contour level value. Must be
                    less than `contours.end`
                usecolormap
                    An alternate to "color". Determines whether or
                    not the contour lines are colored using the
                    trace "colorscale".
                width
                    Sets the width of the contour lines.

        Returns
        -------
        plotly.graph_objs.surface.contours.Z
        """
        return self["z"]

    @z.setter
    def z(self, val):
        self["z"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        x
            :class:`plotly.graph_objects.surface.contours.X`
            instance or dict with compatible properties
        y
            :class:`plotly.graph_objects.surface.contours.Y`
            instance or dict with compatible properties
        z
            :class:`plotly.graph_objects.surface.contours.Z`
            instance or dict with compatible properties
        """

    def __init__(self, arg=None, x=None, y=None, z=None, **kwargs):
        """
        Construct a new Contours object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.surface.Contours`
        x
            :class:`plotly.graph_objects.surface.contours.X`
            instance or dict with compatible properties
        y
            :class:`plotly.graph_objects.surface.contours.Y`
            instance or dict with compatible properties
        z
            :class:`plotly.graph_objects.surface.contours.Z`
            instance or dict with compatible properties

        Returns
        -------
        Contours
        """
        super(Contours, self).__init__("contours")

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
The first argument to the plotly.graph_objs.surface.Contours
constructor must be a dict or
an instance of :class:`plotly.graph_objs.surface.Contours`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("x", None)
        _v = x if x is not None else _v
        if _v is not None:
            self["x"] = _v
        _v = arg.pop("y", None)
        _v = y if y is not None else _v
        if _v is not None:
            self["y"] = _v
        _v = arg.pop("z", None)
        _v = z if z is not None else _v
        if _v is not None:
            self["z"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
