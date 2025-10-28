#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Title(_BaseTraceHierarchyType):
    _parent_path_str = "carpet.baxis"
    _path_str = "carpet.baxis.title"
    _valid_props = {"font", "offset", "text"}

    @property
    def font(self):
        """
        Sets this axis' title font.

        The 'font' property is an instance of Font
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.carpet.baxis.title.Font`
          - A dict of string/value properties that will be passed
            to the Font constructor

        Returns
        -------
        plotly.graph_objs.carpet.baxis.title.Font
        """
        return self["font"]

    @font.setter
    def font(self, val):
        self["font"] = val

    @property
    def offset(self):
        """
        An additional amount by which to offset the title from the tick
        labels, given in pixels.

        The 'offset' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["offset"]

    @offset.setter
    def offset(self, val):
        self["offset"] = val

    @property
    def text(self):
        """
        Sets the title of this axis.

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
            Sets this axis' title font.
        offset
            An additional amount by which to offset the title from
            the tick labels, given in pixels.
        text
            Sets the title of this axis.
        """

    def __init__(self, arg=None, font=None, offset=None, text=None, **kwargs):
        """
        Construct a new Title object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.carpet.baxis.Title`
        font
            Sets this axis' title font.
        offset
            An additional amount by which to offset the title from
            the tick labels, given in pixels.
        text
            Sets the title of this axis.

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
The first argument to the plotly.graph_objs.carpet.baxis.Title
constructor must be a dict or
an instance of :class:`plotly.graph_objs.carpet.baxis.Title`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("font", arg, font)
        self._set_property("offset", arg, offset)
        self._set_property("text", arg, text)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
