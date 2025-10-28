#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Legendgrouptitle(_BaseTraceHierarchyType):
    _parent_path_str = "histogram2d"
    _path_str = "histogram2d.legendgrouptitle"
    _valid_props = {"font", "text"}

    @property
    def font(self):
        """
        Sets this legend group's title font.

        The 'font' property is an instance of Font
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.histogram2d.legendgrouptitle.Font`
          - A dict of string/value properties that will be passed
            to the Font constructor

        Returns
        -------
        plotly.graph_objs.histogram2d.legendgrouptitle.Font
        """
        return self["font"]

    @font.setter
    def font(self, val):
        self["font"] = val

    @property
    def text(self):
        """
        Sets the title of the legend group.

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
            Sets this legend group's title font.
        text
            Sets the title of the legend group.
        """

    def __init__(self, arg=None, font=None, text=None, **kwargs):
        """
        Construct a new Legendgrouptitle object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.histogram2d.Legendgrouptitle`
        font
            Sets this legend group's title font.
        text
            Sets the title of the legend group.

        Returns
        -------
        Legendgrouptitle
        """
        super().__init__("legendgrouptitle")
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
The first argument to the plotly.graph_objs.histogram2d.Legendgrouptitle
constructor must be a dict or
an instance of :class:`plotly.graph_objs.histogram2d.Legendgrouptitle`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("font", arg, font)
        self._set_property("text", arg, text)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
