#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Unselected(_BaseTraceHierarchyType):
    _parent_path_str = "scatterpolargl"
    _path_str = "scatterpolargl.unselected"
    _valid_props = {"marker", "textfont"}

    @property
    def marker(self):
        """
        The 'marker' property is an instance of Marker
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.scatterpolargl.unselected.Marker`
          - A dict of string/value properties that will be passed
            to the Marker constructor

        Returns
        -------
        plotly.graph_objs.scatterpolargl.unselected.Marker
        """
        return self["marker"]

    @marker.setter
    def marker(self, val):
        self["marker"] = val

    @property
    def textfont(self):
        """
        The 'textfont' property is an instance of Textfont
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.scatterpolargl.unselected.Textfont`
          - A dict of string/value properties that will be passed
            to the Textfont constructor

        Returns
        -------
        plotly.graph_objs.scatterpolargl.unselected.Textfont
        """
        return self["textfont"]

    @textfont.setter
    def textfont(self, val):
        self["textfont"] = val

    @property
    def _prop_descriptions(self):
        return """\
        marker
            :class:`plotly.graph_objects.scatterpolargl.unselected.
            Marker` instance or dict with compatible properties
        textfont
            :class:`plotly.graph_objects.scatterpolargl.unselected.
            Textfont` instance or dict with compatible properties
        """

    def __init__(self, arg=None, marker=None, textfont=None, **kwargs):
        """
        Construct a new Unselected object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.scatterpolargl.Unselected`
        marker
            :class:`plotly.graph_objects.scatterpolargl.unselected.
            Marker` instance or dict with compatible properties
        textfont
            :class:`plotly.graph_objects.scatterpolargl.unselected.
            Textfont` instance or dict with compatible properties

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
The first argument to the plotly.graph_objs.scatterpolargl.Unselected
constructor must be a dict or
an instance of :class:`plotly.graph_objs.scatterpolargl.Unselected`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("marker", arg, marker)
        self._set_property("textfont", arg, textfont)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
