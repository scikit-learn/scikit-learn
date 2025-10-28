#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Unifiedhovertitle(_BaseLayoutHierarchyType):
    _parent_path_str = "layout.yaxis"
    _path_str = "layout.yaxis.unifiedhovertitle"
    _valid_props = {"text"}

    @property
    def text(self):
        """
        Template string used for rendering the title that appear on x
        or y unified hover box. Variables are inserted using
        %{variable}, for example "y: %{y}". Numbers are formatted using
        d3-format's syntax %{variable:d3-format}, for example "Price:
        %{y:$.2f}".
        https://github.com/d3/d3-format/tree/v1.4.5#d3-format for
        details on the formatting syntax. Dates are formatted using
        d3-time-format's syntax %{variable|d3-time-format}, for example
        "Day: %{2019-01-01|%A}". https://github.com/d3/d3-time-
        format/tree/v2.2.3#locale_format for details on the date
        formatting syntax.

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
        text
            Template string used for rendering the title that
            appear on x or y unified hover box. Variables are
            inserted using %{variable}, for example "y: %{y}".
            Numbers are formatted using d3-format's syntax
            %{variable:d3-format}, for example "Price: %{y:$.2f}".
            https://github.com/d3/d3-format/tree/v1.4.5#d3-format
            for details on the formatting syntax. Dates are
            formatted using d3-time-format's syntax
            %{variable|d3-time-format}, for example "Day:
            %{2019-01-01|%A}". https://github.com/d3/d3-time-
            format/tree/v2.2.3#locale_format for details on the
            date formatting syntax.
        """

    def __init__(self, arg=None, text=None, **kwargs):
        """
        Construct a new Unifiedhovertitle object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of :class:`plotly.graph_objs.layout.yaxis.U
            nifiedhovertitle`
        text
            Template string used for rendering the title that
            appear on x or y unified hover box. Variables are
            inserted using %{variable}, for example "y: %{y}".
            Numbers are formatted using d3-format's syntax
            %{variable:d3-format}, for example "Price: %{y:$.2f}".
            https://github.com/d3/d3-format/tree/v1.4.5#d3-format
            for details on the formatting syntax. Dates are
            formatted using d3-time-format's syntax
            %{variable|d3-time-format}, for example "Day:
            %{2019-01-01|%A}". https://github.com/d3/d3-time-
            format/tree/v2.2.3#locale_format for details on the
            date formatting syntax.

        Returns
        -------
        Unifiedhovertitle
        """
        super().__init__("unifiedhovertitle")
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
The first argument to the plotly.graph_objs.layout.yaxis.Unifiedhovertitle
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.yaxis.Unifiedhovertitle`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("text", arg, text)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
