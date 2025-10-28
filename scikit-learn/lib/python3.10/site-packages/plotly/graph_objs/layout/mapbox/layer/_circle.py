#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Circle(_BaseLayoutHierarchyType):
    _parent_path_str = "layout.mapbox.layer"
    _path_str = "layout.mapbox.layer.circle"
    _valid_props = {"radius"}

    @property
    def radius(self):
        """
        Sets the circle radius (mapbox.layer.paint.circle-radius). Has
        an effect only when `type` is set to "circle".

        The 'radius' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["radius"]

    @radius.setter
    def radius(self, val):
        self["radius"] = val

    @property
    def _prop_descriptions(self):
        return """\
        radius
            Sets the circle radius (mapbox.layer.paint.circle-
            radius). Has an effect only when `type` is set to
            "circle".
        """

    def __init__(self, arg=None, radius=None, **kwargs):
        """
        Construct a new Circle object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.layout.mapbox.layer.Circle`
        radius
            Sets the circle radius (mapbox.layer.paint.circle-
            radius). Has an effect only when `type` is set to
            "circle".

        Returns
        -------
        Circle
        """
        super().__init__("circle")
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
The first argument to the plotly.graph_objs.layout.mapbox.layer.Circle
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.mapbox.layer.Circle`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("radius", arg, radius)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
