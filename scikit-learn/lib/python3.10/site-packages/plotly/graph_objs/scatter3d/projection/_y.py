#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Y(_BaseTraceHierarchyType):
    _parent_path_str = "scatter3d.projection"
    _path_str = "scatter3d.projection.y"
    _valid_props = {"opacity", "scale", "show"}

    @property
    def opacity(self):
        """
        Sets the projection color.

        The 'opacity' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["opacity"]

    @opacity.setter
    def opacity(self, val):
        self["opacity"] = val

    @property
    def scale(self):
        """
        Sets the scale factor determining the size of the projection
        marker points.

        The 'scale' property is a number and may be specified as:
          - An int or float in the interval [0, 10]

        Returns
        -------
        int|float
        """
        return self["scale"]

    @scale.setter
    def scale(self, val):
        self["scale"] = val

    @property
    def show(self):
        """
        Sets whether or not projections are shown along the y axis.

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

    @property
    def _prop_descriptions(self):
        return """\
        opacity
            Sets the projection color.
        scale
            Sets the scale factor determining the size of the
            projection marker points.
        show
            Sets whether or not projections are shown along the y
            axis.
        """

    def __init__(self, arg=None, opacity=None, scale=None, show=None, **kwargs):
        """
        Construct a new Y object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.scatter3d.projection.Y`
        opacity
            Sets the projection color.
        scale
            Sets the scale factor determining the size of the
            projection marker points.
        show
            Sets whether or not projections are shown along the y
            axis.

        Returns
        -------
        Y
        """
        super().__init__("y")
        if "_parent" in kwargs:
            self._parent = kwargs["_parent"]
            return

        if arg is None:
            arg = {}
        elif isinstance(arg, self.__class__):
            arg = arg.to_plotly_json()
        elif isinstance(arg, dict):
            arg = _copy.copy(arg)
        else:
            raise ValueError("""\
The first argument to the plotly.graph_objs.scatter3d.projection.Y
constructor must be a dict or
an instance of :class:`plotly.graph_objs.scatter3d.projection.Y`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("opacity", arg, opacity)
        self._set_property("scale", arg, scale)
        self._set_property("show", arg, show)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
