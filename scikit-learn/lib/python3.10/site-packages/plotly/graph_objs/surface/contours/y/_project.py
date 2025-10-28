#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Project(_BaseTraceHierarchyType):
    _parent_path_str = "surface.contours.y"
    _path_str = "surface.contours.y.project"
    _valid_props = {"x", "y", "z"}

    @property
    def x(self):
        """
        Determines whether or not these contour lines are projected on
        the x plane. If `highlight` is set to True (the default), the
        projected lines are shown on hover. If `show` is set to True,
        the projected lines are shown in permanence.

        The 'x' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["x"]

    @x.setter
    def x(self, val):
        self["x"] = val

    @property
    def y(self):
        """
        Determines whether or not these contour lines are projected on
        the y plane. If `highlight` is set to True (the default), the
        projected lines are shown on hover. If `show` is set to True,
        the projected lines are shown in permanence.

        The 'y' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["y"]

    @y.setter
    def y(self, val):
        self["y"] = val

    @property
    def z(self):
        """
        Determines whether or not these contour lines are projected on
        the z plane. If `highlight` is set to True (the default), the
        projected lines are shown on hover. If `show` is set to True,
        the projected lines are shown in permanence.

        The 'z' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["z"]

    @z.setter
    def z(self, val):
        self["z"] = val

    @property
    def _prop_descriptions(self):
        return """\
        x
            Determines whether or not these contour lines are
            projected on the x plane. If `highlight` is set to True
            (the default), the projected lines are shown on hover.
            If `show` is set to True, the projected lines are shown
            in permanence.
        y
            Determines whether or not these contour lines are
            projected on the y plane. If `highlight` is set to True
            (the default), the projected lines are shown on hover.
            If `show` is set to True, the projected lines are shown
            in permanence.
        z
            Determines whether or not these contour lines are
            projected on the z plane. If `highlight` is set to True
            (the default), the projected lines are shown on hover.
            If `show` is set to True, the projected lines are shown
            in permanence.
        """

    def __init__(self, arg=None, x=None, y=None, z=None, **kwargs):
        """
        Construct a new Project object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.surface.contours.y.Project`
        x
            Determines whether or not these contour lines are
            projected on the x plane. If `highlight` is set to True
            (the default), the projected lines are shown on hover.
            If `show` is set to True, the projected lines are shown
            in permanence.
        y
            Determines whether or not these contour lines are
            projected on the y plane. If `highlight` is set to True
            (the default), the projected lines are shown on hover.
            If `show` is set to True, the projected lines are shown
            in permanence.
        z
            Determines whether or not these contour lines are
            projected on the z plane. If `highlight` is set to True
            (the default), the projected lines are shown on hover.
            If `show` is set to True, the projected lines are shown
            in permanence.

        Returns
        -------
        Project
        """
        super().__init__("project")
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
The first argument to the plotly.graph_objs.surface.contours.y.Project
constructor must be a dict or
an instance of :class:`plotly.graph_objs.surface.contours.y.Project`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("x", arg, x)
        self._set_property("y", arg, y)
        self._set_property("z", arg, z)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
