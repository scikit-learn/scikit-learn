#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Fill(_BaseLayoutHierarchyType):
    _parent_path_str = "layout.map.layer"
    _path_str = "layout.map.layer.fill"
    _valid_props = {"outlinecolor"}

    @property
    def outlinecolor(self):
        """
        Sets the fill outline color (map.layer.paint.fill-outline-
        color). Has an effect only when `type` is set to "fill".

        The 'outlinecolor' property is a color and may be specified as:
          - A hex string (e.g. '#ff0000')
          - An rgb/rgba string (e.g. 'rgb(255,0,0)')
          - An hsl/hsla string (e.g. 'hsl(0,100%,50%)')
          - An hsv/hsva string (e.g. 'hsv(0,100%,100%)')
          - A named CSS color: see https://plotly.com/python/css-colors/ for a list

        Returns
        -------
        str
        """
        return self["outlinecolor"]

    @outlinecolor.setter
    def outlinecolor(self, val):
        self["outlinecolor"] = val

    @property
    def _prop_descriptions(self):
        return """\
        outlinecolor
            Sets the fill outline color (map.layer.paint.fill-
            outline-color). Has an effect only when `type` is set
            to "fill".
        """

    def __init__(self, arg=None, outlinecolor=None, **kwargs):
        """
        Construct a new Fill object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.layout.map.layer.Fill`
        outlinecolor
            Sets the fill outline color (map.layer.paint.fill-
            outline-color). Has an effect only when `type` is set
            to "fill".

        Returns
        -------
        Fill
        """
        super().__init__("fill")
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
The first argument to the plotly.graph_objs.layout.map.layer.Fill
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.map.layer.Fill`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("outlinecolor", arg, outlinecolor)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
