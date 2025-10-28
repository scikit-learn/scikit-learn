#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Unselected(_BaseTraceHierarchyType):
    _parent_path_str = "parcoords"
    _path_str = "parcoords.unselected"
    _valid_props = {"line"}

    @property
    def line(self):
        """
        The 'line' property is an instance of Line
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.parcoords.unselected.Line`
          - A dict of string/value properties that will be passed
            to the Line constructor

        Returns
        -------
        plotly.graph_objs.parcoords.unselected.Line
        """
        return self["line"]

    @line.setter
    def line(self, val):
        self["line"] = val

    @property
    def _prop_descriptions(self):
        return """\
        line
            :class:`plotly.graph_objects.parcoords.unselected.Line`
            instance or dict with compatible properties
        """

    def __init__(self, arg=None, line=None, **kwargs):
        """
        Construct a new Unselected object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.parcoords.Unselected`
        line
            :class:`plotly.graph_objects.parcoords.unselected.Line`
            instance or dict with compatible properties

        Returns
        -------
        Unselected
        """
        super().__init__("unselected")
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
The first argument to the plotly.graph_objs.parcoords.Unselected
constructor must be a dict or
an instance of :class:`plotly.graph_objs.parcoords.Unselected`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("line", arg, line)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
