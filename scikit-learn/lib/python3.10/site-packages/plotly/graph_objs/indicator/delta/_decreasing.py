#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Decreasing(_BaseTraceHierarchyType):
    _parent_path_str = "indicator.delta"
    _path_str = "indicator.delta.decreasing"
    _valid_props = {"color", "symbol"}

    @property
    def color(self):
        """
        Sets the color for increasing value.

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
    def symbol(self):
        """
        Sets the symbol to display for increasing value

        The 'symbol' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["symbol"]

    @symbol.setter
    def symbol(self, val):
        self["symbol"] = val

    @property
    def _prop_descriptions(self):
        return """\
        color
            Sets the color for increasing value.
        symbol
            Sets the symbol to display for increasing value
        """

    def __init__(self, arg=None, color=None, symbol=None, **kwargs):
        """
        Construct a new Decreasing object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.indicator.delta.Decreasing`
        color
            Sets the color for increasing value.
        symbol
            Sets the symbol to display for increasing value

        Returns
        -------
        Decreasing
        """
        super().__init__("decreasing")
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
The first argument to the plotly.graph_objs.indicator.delta.Decreasing
constructor must be a dict or
an instance of :class:`plotly.graph_objs.indicator.delta.Decreasing`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("color", arg, color)
        self._set_property("symbol", arg, symbol)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
