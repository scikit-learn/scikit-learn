#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Title(_BaseLayoutHierarchyType):
    _parent_path_str = "layout.legend"
    _path_str = "layout.legend.title"
    _valid_props = {"font", "side", "text"}

    @property
    def font(self):
        """
        Sets this legend's title font. Defaults to `legend.font` with
        its size increased about 20%.

        The 'font' property is an instance of Font
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.legend.title.Font`
          - A dict of string/value properties that will be passed
            to the Font constructor

        Returns
        -------
        plotly.graph_objs.layout.legend.title.Font
        """
        return self["font"]

    @font.setter
    def font(self, val):
        self["font"] = val

    @property
    def side(self):
        """
        Determines the location of legend's title with respect to the
        legend items. Defaulted to "top" with `orientation` is "h".
        Defaulted to "left" with `orientation` is "v". The *top left*
        options could be used to expand top center and top right are
        for horizontal alignment legend area in both x and y sides.

        The 'side' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['top', 'left', 'top left', 'top center', 'top right']

        Returns
        -------
        Any
        """
        return self["side"]

    @side.setter
    def side(self, val):
        self["side"] = val

    @property
    def text(self):
        """
        Sets the title of the legend.

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
    def _prop_descriptions(self):
        return """\
        font
            Sets this legend's title font. Defaults to
            `legend.font` with its size increased about 20%.
        side
            Determines the location of legend's title with respect
            to the legend items. Defaulted to "top" with
            `orientation` is "h". Defaulted to "left" with
            `orientation` is "v". The *top left* options could be
            used to expand top center and top right are for
            horizontal alignment legend area in both x and y sides.
        text
            Sets the title of the legend.
        """

    def __init__(self, arg=None, font=None, side=None, text=None, **kwargs):
        """
        Construct a new Title object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.layout.legend.Title`
        font
            Sets this legend's title font. Defaults to
            `legend.font` with its size increased about 20%.
        side
            Determines the location of legend's title with respect
            to the legend items. Defaulted to "top" with
            `orientation` is "h". Defaulted to "left" with
            `orientation` is "v". The *top left* options could be
            used to expand top center and top right are for
            horizontal alignment legend area in both x and y sides.
        text
            Sets the title of the legend.

        Returns
        -------
        Title
        """
        super().__init__("title")
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
The first argument to the plotly.graph_objs.layout.legend.Title
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.legend.Title`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("font", arg, font)
        self._set_property("side", arg, side)
        self._set_property("text", arg, text)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
