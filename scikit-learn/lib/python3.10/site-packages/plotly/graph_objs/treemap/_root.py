#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Root(_BaseTraceHierarchyType):
    _parent_path_str = "treemap"
    _path_str = "treemap.root"
    _valid_props = {"color"}

    @property
    def color(self):
        """
        sets the color of the root node for a sunburst/treemap/icicle
        trace. this has no effect when a colorscale is used to set the
        markers.

        The 'color' property is a color and may be specified as:
          - A hex string (e.g. '#ff0000')
          - An rgb/rgba string (e.g. 'rgb(255,0,0)')
          - An hsl/hsla string (e.g. 'hsl(0,100%,50%)')
          - An hsv/hsva string (e.g. 'hsv(0,100%,100%)')
          - A named CSS color: see https://plotly.com/python/css-colors/ for a list

        Returns
        -------
        str
        """
        return self["color"]

    @color.setter
    def color(self, val):
        self["color"] = val

    @property
    def _prop_descriptions(self):
        return """\
        color
            sets the color of the root node for a
            sunburst/treemap/icicle trace. this has no effect when
            a colorscale is used to set the markers.
        """

    def __init__(self, arg=None, color=None, **kwargs):
        """
        Construct a new Root object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of :class:`plotly.graph_objs.treemap.Root`
        color
            sets the color of the root node for a
            sunburst/treemap/icicle trace. this has no effect when
            a colorscale is used to set the markers.

        Returns
        -------
        Root
        """
        super().__init__("root")
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
The first argument to the plotly.graph_objs.treemap.Root
constructor must be a dict or
an instance of :class:`plotly.graph_objs.treemap.Root`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("color", arg, color)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
