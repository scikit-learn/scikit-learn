#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Aspectratio(_BaseLayoutHierarchyType):
    _parent_path_str = "layout.scene"
    _path_str = "layout.scene.aspectratio"
    _valid_props = {"x", "y", "z"}

    @property
    def x(self):
        """
        The 'x' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["x"]

    @x.setter
    def x(self, val):
        self["x"] = val

    @property
    def y(self):
        """
        The 'y' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["y"]

    @y.setter
    def y(self, val):
        self["y"] = val

    @property
    def z(self):
        """
        The 'z' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["z"]

    @z.setter
    def z(self, val):
        self["z"] = val

    @property
    def _prop_descriptions(self):
        return """\
        x

        y

        z

        """

    def __init__(self, arg=None, x=None, y=None, z=None, **kwargs):
        """
        Construct a new Aspectratio object

        Sets this scene's axis aspectratio.

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.layout.scene.Aspectratio`
        x

        y

        z


        Returns
        -------
        Aspectratio
        """
        super().__init__("aspectratio")
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
The first argument to the plotly.graph_objs.layout.scene.Aspectratio
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.scene.Aspectratio`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("x", arg, x)
        self._set_property("y", arg, y)
        self._set_property("z", arg, z)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
