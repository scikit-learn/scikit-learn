#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Title(_BaseTraceHierarchyType):
    _parent_path_str = "scatter.marker.colorbar"
    _path_str = "scatter.marker.colorbar.title"
    _valid_props = {"font", "side", "text"}

    @property
    def font(self):
        """
        Sets this color bar's title font.

        The 'font' property is an instance of Font
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.scatter.marker.colorbar.title.Font`
          - A dict of string/value properties that will be passed
            to the Font constructor

        Returns
        -------
        plotly.graph_objs.scatter.marker.colorbar.title.Font
        """
        return self["font"]

    @font.setter
    def font(self, val):
        self["font"] = val

    @property
    def side(self):
        """
        Determines the location of color bar's title with respect to
        the color bar. Defaults to "top" when `orientation` if "v" and
        defaults to "right" when `orientation` if "h".

        The 'side' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['right', 'top', 'bottom']

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
        Sets the title of the color bar.

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
            Sets this color bar's title font.
        side
            Determines the location of color bar's title with
            respect to the color bar. Defaults to "top" when
            `orientation` if "v" and  defaults to "right" when
            `orientation` if "h".
        text
            Sets the title of the color bar.
        """

    def __init__(self, arg=None, font=None, side=None, text=None, **kwargs):
        """
        Construct a new Title object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of :class:`plotly.graph_objs.scatter.marker
            .colorbar.Title`
        font
            Sets this color bar's title font.
        side
            Determines the location of color bar's title with
            respect to the color bar. Defaults to "top" when
            `orientation` if "v" and  defaults to "right" when
            `orientation` if "h".
        text
            Sets the title of the color bar.

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
The first argument to the plotly.graph_objs.scatter.marker.colorbar.Title
constructor must be a dict or
an instance of :class:`plotly.graph_objs.scatter.marker.colorbar.Title`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("font", arg, font)
        self._set_property("side", arg, side)
        self._set_property("text", arg, text)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
