from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Slices(_BaseTraceHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "isosurface"
    _path_str = "isosurface.slices"
    _valid_props = {"x", "y", "z"}

    # x
    # -
    @property
    def x(self):
        """
        The 'x' property is an instance of X
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.isosurface.slices.X`
          - A dict of string/value properties that will be passed
            to the X constructor

            Supported dict properties:

                fill
                    Sets the fill ratio of the `slices`. The
                    default fill value of the `slices` is 1 meaning
                    that they are entirely shaded. On the other
                    hand Applying a `fill` ratio less than one
                    would allow the creation of openings parallel
                    to the edges.
                locations
                    Specifies the location(s) of slices on the
                    axis. When not specified slices would be
                    created for all points of the axis x except
                    start and end.
                locationssrc
                    Sets the source reference on Chart Studio Cloud
                    for `locations`.
                show
                    Determines whether or not slice planes about
                    the x dimension are drawn.

        Returns
        -------
        plotly.graph_objs.isosurface.slices.X
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
          - An instance of :class:`plotly.graph_objs.isosurface.slices.Y`
          - A dict of string/value properties that will be passed
            to the Y constructor

            Supported dict properties:

                fill
                    Sets the fill ratio of the `slices`. The
                    default fill value of the `slices` is 1 meaning
                    that they are entirely shaded. On the other
                    hand Applying a `fill` ratio less than one
                    would allow the creation of openings parallel
                    to the edges.
                locations
                    Specifies the location(s) of slices on the
                    axis. When not specified slices would be
                    created for all points of the axis y except
                    start and end.
                locationssrc
                    Sets the source reference on Chart Studio Cloud
                    for `locations`.
                show
                    Determines whether or not slice planes about
                    the y dimension are drawn.

        Returns
        -------
        plotly.graph_objs.isosurface.slices.Y
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
          - An instance of :class:`plotly.graph_objs.isosurface.slices.Z`
          - A dict of string/value properties that will be passed
            to the Z constructor

            Supported dict properties:

                fill
                    Sets the fill ratio of the `slices`. The
                    default fill value of the `slices` is 1 meaning
                    that they are entirely shaded. On the other
                    hand Applying a `fill` ratio less than one
                    would allow the creation of openings parallel
                    to the edges.
                locations
                    Specifies the location(s) of slices on the
                    axis. When not specified slices would be
                    created for all points of the axis z except
                    start and end.
                locationssrc
                    Sets the source reference on Chart Studio Cloud
                    for `locations`.
                show
                    Determines whether or not slice planes about
                    the z dimension are drawn.

        Returns
        -------
        plotly.graph_objs.isosurface.slices.Z
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
            :class:`plotly.graph_objects.isosurface.slices.X`
            instance or dict with compatible properties
        y
            :class:`plotly.graph_objects.isosurface.slices.Y`
            instance or dict with compatible properties
        z
            :class:`plotly.graph_objects.isosurface.slices.Z`
            instance or dict with compatible properties
        """

    def __init__(self, arg=None, x=None, y=None, z=None, **kwargs):
        """
        Construct a new Slices object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.isosurface.Slices`
        x
            :class:`plotly.graph_objects.isosurface.slices.X`
            instance or dict with compatible properties
        y
            :class:`plotly.graph_objects.isosurface.slices.Y`
            instance or dict with compatible properties
        z
            :class:`plotly.graph_objects.isosurface.slices.Z`
            instance or dict with compatible properties

        Returns
        -------
        Slices
        """
        super(Slices, self).__init__("slices")

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
The first argument to the plotly.graph_objs.isosurface.Slices
constructor must be a dict or
an instance of :class:`plotly.graph_objs.isosurface.Slices`"""
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
