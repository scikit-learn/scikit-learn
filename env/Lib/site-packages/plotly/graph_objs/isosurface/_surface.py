from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Surface(_BaseTraceHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "isosurface"
    _path_str = "isosurface.surface"
    _valid_props = {"count", "fill", "pattern", "show"}

    # count
    # -----
    @property
    def count(self):
        """
        Sets the number of iso-surfaces between minimum and maximum
        iso-values. By default this value is 2 meaning that only
        minimum and maximum surfaces would be drawn.

        The 'count' property is a integer and may be specified as:
          - An int (or float that will be cast to an int)
            in the interval [1, 9223372036854775807]

        Returns
        -------
        int
        """
        return self["count"]

    @count.setter
    def count(self, val):
        self["count"] = val

    # fill
    # ----
    @property
    def fill(self):
        """
        Sets the fill ratio of the iso-surface. The default fill value
        of the surface is 1 meaning that they are entirely shaded. On
        the other hand Applying a `fill` ratio less than one would
        allow the creation of openings parallel to the edges.

        The 'fill' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["fill"]

    @fill.setter
    def fill(self, val):
        self["fill"] = val

    # pattern
    # -------
    @property
    def pattern(self):
        """
        Sets the surface pattern of the iso-surface 3-D sections. The
        default pattern of the surface is `all` meaning that the rest
        of surface elements would be shaded. The check options (either
        1 or 2) could be used to draw half of the squares on the
        surface. Using various combinations of capital `A`, `B`, `C`,
        `D` and `E` may also be used to reduce the number of triangles
        on the iso-surfaces and creating other patterns of interest.

        The 'pattern' property is a flaglist and may be specified
        as a string containing:
          - Any combination of ['A', 'B', 'C', 'D', 'E'] joined with '+' characters
            (e.g. 'A+B')
            OR exactly one of ['all', 'odd', 'even'] (e.g. 'even')

        Returns
        -------
        Any
        """
        return self["pattern"]

    @pattern.setter
    def pattern(self, val):
        self["pattern"] = val

    # show
    # ----
    @property
    def show(self):
        """
        Hides/displays surfaces between minimum and maximum iso-values.

        The 'show' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["show"]

    @show.setter
    def show(self, val):
        self["show"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        count
            Sets the number of iso-surfaces between minimum and
            maximum iso-values. By default this value is 2 meaning
            that only minimum and maximum surfaces would be drawn.
        fill
            Sets the fill ratio of the iso-surface. The default
            fill value of the surface is 1 meaning that they are
            entirely shaded. On the other hand Applying a `fill`
            ratio less than one would allow the creation of
            openings parallel to the edges.
        pattern
            Sets the surface pattern of the iso-surface 3-D
            sections. The default pattern of the surface is `all`
            meaning that the rest of surface elements would be
            shaded. The check options (either 1 or 2) could be used
            to draw half of the squares on the surface. Using
            various combinations of capital `A`, `B`, `C`, `D` and
            `E` may also be used to reduce the number of triangles
            on the iso-surfaces and creating other patterns of
            interest.
        show
            Hides/displays surfaces between minimum and maximum
            iso-values.
        """

    def __init__(
        self, arg=None, count=None, fill=None, pattern=None, show=None, **kwargs
    ):
        """
        Construct a new Surface object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.isosurface.Surface`
        count
            Sets the number of iso-surfaces between minimum and
            maximum iso-values. By default this value is 2 meaning
            that only minimum and maximum surfaces would be drawn.
        fill
            Sets the fill ratio of the iso-surface. The default
            fill value of the surface is 1 meaning that they are
            entirely shaded. On the other hand Applying a `fill`
            ratio less than one would allow the creation of
            openings parallel to the edges.
        pattern
            Sets the surface pattern of the iso-surface 3-D
            sections. The default pattern of the surface is `all`
            meaning that the rest of surface elements would be
            shaded. The check options (either 1 or 2) could be used
            to draw half of the squares on the surface. Using
            various combinations of capital `A`, `B`, `C`, `D` and
            `E` may also be used to reduce the number of triangles
            on the iso-surfaces and creating other patterns of
            interest.
        show
            Hides/displays surfaces between minimum and maximum
            iso-values.

        Returns
        -------
        Surface
        """
        super(Surface, self).__init__("surface")

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
The first argument to the plotly.graph_objs.isosurface.Surface
constructor must be a dict or
an instance of :class:`plotly.graph_objs.isosurface.Surface`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("count", None)
        _v = count if count is not None else _v
        if _v is not None:
            self["count"] = _v
        _v = arg.pop("fill", None)
        _v = fill if fill is not None else _v
        if _v is not None:
            self["fill"] = _v
        _v = arg.pop("pattern", None)
        _v = pattern if pattern is not None else _v
        if _v is not None:
            self["pattern"] = _v
        _v = arg.pop("show", None)
        _v = show if show is not None else _v
        if _v is not None:
            self["show"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
