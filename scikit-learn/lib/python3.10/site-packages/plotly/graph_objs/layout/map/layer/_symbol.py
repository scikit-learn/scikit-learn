#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Symbol(_BaseLayoutHierarchyType):
    _parent_path_str = "layout.map.layer"
    _path_str = "layout.map.layer.symbol"
    _valid_props = {"icon", "iconsize", "placement", "text", "textfont", "textposition"}

    @property
    def icon(self):
        """
        Sets the symbol icon image (map.layer.layout.icon-image). Full
        list: https://www.mapbox.com/maki-icons/

        The 'icon' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["icon"]

    @icon.setter
    def icon(self, val):
        self["icon"] = val

    @property
    def iconsize(self):
        """
        Sets the symbol icon size (map.layer.layout.icon-size). Has an
        effect only when `type` is set to "symbol".

        The 'iconsize' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["iconsize"]

    @iconsize.setter
    def iconsize(self, val):
        self["iconsize"] = val

    @property
    def placement(self):
        """
        Sets the symbol and/or text placement (map.layer.layout.symbol-
        placement). If `placement` is "point", the label is placed
        where the geometry is located If `placement` is "line", the
        label is placed along the line of the geometry If `placement`
        is "line-center", the label is placed on the center of the
        geometry

        The 'placement' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['point', 'line', 'line-center']

        Returns
        -------
        Any
        """
        return self["placement"]

    @placement.setter
    def placement(self, val):
        self["placement"] = val

    @property
    def text(self):
        """
        Sets the symbol text (map.layer.layout.text-field).

        The 'text' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["text"]

    @text.setter
    def text(self, val):
        self["text"] = val

    @property
    def textfont(self):
        """
        Sets the icon text font (color=map.layer.paint.text-color,
        size=map.layer.layout.text-size). Has an effect only when
        `type` is set to "symbol".

        The 'textfont' property is an instance of Textfont
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.map.layer.symbol.Textfont`
          - A dict of string/value properties that will be passed
            to the Textfont constructor

        Returns
        -------
        plotly.graph_objs.layout.map.layer.symbol.Textfont
        """
        return self["textfont"]

    @textfont.setter
    def textfont(self, val):
        self["textfont"] = val

    @property
    def textposition(self):
        """
        Sets the positions of the `text` elements with respects to the
        (x,y) coordinates.

        The 'textposition' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['top left', 'top center', 'top right', 'middle left',
                'middle center', 'middle right', 'bottom left', 'bottom
                center', 'bottom right']

        Returns
        -------
        Any
        """
        return self["textposition"]

    @textposition.setter
    def textposition(self, val):
        self["textposition"] = val

    @property
    def _prop_descriptions(self):
        return """\
        icon
            Sets the symbol icon image (map.layer.layout.icon-
            image). Full list: https://www.mapbox.com/maki-icons/
        iconsize
            Sets the symbol icon size (map.layer.layout.icon-size).
            Has an effect only when `type` is set to "symbol".
        placement
            Sets the symbol and/or text placement
            (map.layer.layout.symbol-placement). If `placement` is
            "point", the label is placed where the geometry is
            located If `placement` is "line", the label is placed
            along the line of the geometry If `placement` is "line-
            center", the label is placed on the center of the
            geometry
        text
            Sets the symbol text (map.layer.layout.text-field).
        textfont
            Sets the icon text font (color=map.layer.paint.text-
            color, size=map.layer.layout.text-size). Has an effect
            only when `type` is set to "symbol".
        textposition
            Sets the positions of the `text` elements with respects
            to the (x,y) coordinates.
        """

    def __init__(
        self,
        arg=None,
        icon=None,
        iconsize=None,
        placement=None,
        text=None,
        textfont=None,
        textposition=None,
        **kwargs,
    ):
        """
        Construct a new Symbol object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.layout.map.layer.Symbol`
        icon
            Sets the symbol icon image (map.layer.layout.icon-
            image). Full list: https://www.mapbox.com/maki-icons/
        iconsize
            Sets the symbol icon size (map.layer.layout.icon-size).
            Has an effect only when `type` is set to "symbol".
        placement
            Sets the symbol and/or text placement
            (map.layer.layout.symbol-placement). If `placement` is
            "point", the label is placed where the geometry is
            located If `placement` is "line", the label is placed
            along the line of the geometry If `placement` is "line-
            center", the label is placed on the center of the
            geometry
        text
            Sets the symbol text (map.layer.layout.text-field).
        textfont
            Sets the icon text font (color=map.layer.paint.text-
            color, size=map.layer.layout.text-size). Has an effect
            only when `type` is set to "symbol".
        textposition
            Sets the positions of the `text` elements with respects
            to the (x,y) coordinates.

        Returns
        -------
        Symbol
        """
        super().__init__("symbol")
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
The first argument to the plotly.graph_objs.layout.map.layer.Symbol
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.map.layer.Symbol`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("icon", arg, icon)
        self._set_property("iconsize", arg, iconsize)
        self._set_property("placement", arg, placement)
        self._set_property("text", arg, text)
        self._set_property("textfont", arg, textfont)
        self._set_property("textposition", arg, textposition)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
