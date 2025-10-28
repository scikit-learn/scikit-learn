#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Caps(_BaseTraceHierarchyType):
    _parent_path_str = "isosurface"
    _path_str = "isosurface.caps"
    _valid_props = {"x", "y", "z"}

    @property
    def x(self):
        """
        The 'x' property is an instance of X
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.isosurface.caps.X`
          - A dict of string/value properties that will be passed
            to the X constructor

        Returns
        -------
        plotly.graph_objs.isosurface.caps.X
        """
        return self["x"]

    @x.setter
    def x(self, val):
        self["x"] = val

    @property
    def y(self):
        """
        The 'y' property is an instance of Y
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.isosurface.caps.Y`
          - A dict of string/value properties that will be passed
            to the Y constructor

        Returns
        -------
        plotly.graph_objs.isosurface.caps.Y
        """
        return self["y"]

    @y.setter
    def y(self, val):
        self["y"] = val

    @property
    def z(self):
        """
        The 'z' property is an instance of Z
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.isosurface.caps.Z`
          - A dict of string/value properties that will be passed
            to the Z constructor

        Returns
        -------
        plotly.graph_objs.isosurface.caps.Z
        """
        return self["z"]

    @z.setter
    def z(self, val):
        self["z"] = val

    @property
    def _prop_descriptions(self):
        return """\
        x
            :class:`plotly.graph_objects.isosurface.caps.X`
            instance or dict with compatible properties
        y
            :class:`plotly.graph_objects.isosurface.caps.Y`
            instance or dict with compatible properties
        z
            :class:`plotly.graph_objects.isosurface.caps.Z`
            instance or dict with compatible properties
        """

    def __init__(self, arg=None, x=None, y=None, z=None, **kwargs):
        """
        Construct a new Caps object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.isosurface.Caps`
        x
            :class:`plotly.graph_objects.isosurface.caps.X`
            instance or dict with compatible properties
        y
            :class:`plotly.graph_objects.isosurface.caps.Y`
            instance or dict with compatible properties
        z
            :class:`plotly.graph_objects.isosurface.caps.Z`
            instance or dict with compatible properties

        Returns
        -------
        Caps
        """
        super().__init__("caps")
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
The first argument to the plotly.graph_objs.isosurface.Caps
constructor must be a dict or
an instance of :class:`plotly.graph_objs.isosurface.Caps`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("x", arg, x)
        self._set_property("y", arg, y)
        self._set_property("z", arg, z)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
