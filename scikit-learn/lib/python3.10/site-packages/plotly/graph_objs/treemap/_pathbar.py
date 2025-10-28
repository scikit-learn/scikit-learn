#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Pathbar(_BaseTraceHierarchyType):
    _parent_path_str = "treemap"
    _path_str = "treemap.pathbar"
    _valid_props = {"edgeshape", "side", "textfont", "thickness", "visible"}

    @property
    def edgeshape(self):
        """
        Determines which shape is used for edges between `barpath`
        labels.

        The 'edgeshape' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['>', '<', '|', '/', '\\']

        Returns
        -------
        Any
        """
        return self["edgeshape"]

    @edgeshape.setter
    def edgeshape(self, val):
        self["edgeshape"] = val

    @property
    def side(self):
        """
        Determines on which side of the the treemap the `pathbar`
        should be presented.

        The 'side' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['top', 'bottom']

        Returns
        -------
        Any
        """
        return self["side"]

    @side.setter
    def side(self, val):
        self["side"] = val

    @property
    def textfont(self):
        """
        Sets the font used inside `pathbar`.

        The 'textfont' property is an instance of Textfont
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.treemap.pathbar.Textfont`
          - A dict of string/value properties that will be passed
            to the Textfont constructor

        Returns
        -------
        plotly.graph_objs.treemap.pathbar.Textfont
        """
        return self["textfont"]

    @textfont.setter
    def textfont(self, val):
        self["textfont"] = val

    @property
    def thickness(self):
        """
        Sets the thickness of `pathbar` (in px). If not specified the
        `pathbar.textfont.size` is used with 3 pixles extra padding on
        each side.

        The 'thickness' property is a number and may be specified as:
          - An int or float in the interval [12, inf]

        Returns
        -------
        int|float
        """
        return self["thickness"]

    @thickness.setter
    def thickness(self, val):
        self["thickness"] = val

    @property
    def visible(self):
        """
        Determines if the path bar is drawn i.e. outside the trace
        `domain` and with one pixel gap.

        The 'visible' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["visible"]

    @visible.setter
    def visible(self, val):
        self["visible"] = val

    @property
    def _prop_descriptions(self):
        return """\
        edgeshape
            Determines which shape is used for edges between
            `barpath` labels.
        side
            Determines on which side of the the treemap the
            `pathbar` should be presented.
        textfont
            Sets the font used inside `pathbar`.
        thickness
            Sets the thickness of `pathbar` (in px). If not
            specified the `pathbar.textfont.size` is used with 3
            pixles extra padding on each side.
        visible
            Determines if the path bar is drawn i.e. outside the
            trace `domain` and with one pixel gap.
        """

    def __init__(
        self,
        arg=None,
        edgeshape=None,
        side=None,
        textfont=None,
        thickness=None,
        visible=None,
        **kwargs,
    ):
        """
        Construct a new Pathbar object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.treemap.Pathbar`
        edgeshape
            Determines which shape is used for edges between
            `barpath` labels.
        side
            Determines on which side of the the treemap the
            `pathbar` should be presented.
        textfont
            Sets the font used inside `pathbar`.
        thickness
            Sets the thickness of `pathbar` (in px). If not
            specified the `pathbar.textfont.size` is used with 3
            pixles extra padding on each side.
        visible
            Determines if the path bar is drawn i.e. outside the
            trace `domain` and with one pixel gap.

        Returns
        -------
        Pathbar
        """
        super().__init__("pathbar")
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
The first argument to the plotly.graph_objs.treemap.Pathbar
constructor must be a dict or
an instance of :class:`plotly.graph_objs.treemap.Pathbar`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("edgeshape", arg, edgeshape)
        self._set_property("side", arg, side)
        self._set_property("textfont", arg, textfont)
        self._set_property("thickness", arg, thickness)
        self._set_property("visible", arg, visible)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
