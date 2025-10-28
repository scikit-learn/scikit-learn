#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Newselection(_BaseLayoutHierarchyType):
    _parent_path_str = "layout"
    _path_str = "layout.newselection"
    _valid_props = {"line", "mode"}

    @property
    def line(self):
        """
        The 'line' property is an instance of Line
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.newselection.Line`
          - A dict of string/value properties that will be passed
            to the Line constructor

        Returns
        -------
        plotly.graph_objs.layout.newselection.Line
        """
        return self["line"]

    @line.setter
    def line(self, val):
        self["line"] = val

    @property
    def mode(self):
        """
        Describes how a new selection is created. If `immediate`, a new
        selection is created after first mouse up. If `gradual`, a new
        selection is not created after first mouse. By adding to and
        subtracting from the initial selection, this option allows
        declaring extra outlines of the selection.

        The 'mode' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['immediate', 'gradual']

        Returns
        -------
        Any
        """
        return self["mode"]

    @mode.setter
    def mode(self, val):
        self["mode"] = val

    @property
    def _prop_descriptions(self):
        return """\
        line
            :class:`plotly.graph_objects.layout.newselection.Line`
            instance or dict with compatible properties
        mode
            Describes how a new selection is created. If
            `immediate`, a new selection is created after first
            mouse up. If `gradual`, a new selection is not created
            after first mouse. By adding to and subtracting from
            the initial selection, this option allows declaring
            extra outlines of the selection.
        """

    def __init__(self, arg=None, line=None, mode=None, **kwargs):
        """
        Construct a new Newselection object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.layout.Newselection`
        line
            :class:`plotly.graph_objects.layout.newselection.Line`
            instance or dict with compatible properties
        mode
            Describes how a new selection is created. If
            `immediate`, a new selection is created after first
            mouse up. If `gradual`, a new selection is not created
            after first mouse. By adding to and subtracting from
            the initial selection, this option allows declaring
            extra outlines of the selection.

        Returns
        -------
        Newselection
        """
        super().__init__("newselection")
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
The first argument to the plotly.graph_objs.layout.Newselection
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.Newselection`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("line", arg, line)
        self._set_property("mode", arg, mode)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
