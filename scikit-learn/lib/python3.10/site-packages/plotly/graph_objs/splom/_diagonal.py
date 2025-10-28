#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Diagonal(_BaseTraceHierarchyType):
    _parent_path_str = "splom"
    _path_str = "splom.diagonal"
    _valid_props = {"visible"}

    @property
    def visible(self):
        """
        Determines whether or not subplots on the diagonal are
        displayed.

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
        visible
            Determines whether or not subplots on the diagonal are
            displayed.
        """

    def __init__(self, arg=None, visible=None, **kwargs):
        """
        Construct a new Diagonal object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.splom.Diagonal`
        visible
            Determines whether or not subplots on the diagonal are
            displayed.

        Returns
        -------
        Diagonal
        """
        super().__init__("diagonal")
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
The first argument to the plotly.graph_objs.splom.Diagonal
constructor must be a dict or
an instance of :class:`plotly.graph_objs.splom.Diagonal`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("visible", arg, visible)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
