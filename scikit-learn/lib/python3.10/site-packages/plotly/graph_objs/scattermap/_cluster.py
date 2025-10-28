#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Cluster(_BaseTraceHierarchyType):
    _parent_path_str = "scattermap"
    _path_str = "scattermap.cluster"
    _valid_props = {
        "color",
        "colorsrc",
        "enabled",
        "maxzoom",
        "opacity",
        "opacitysrc",
        "size",
        "sizesrc",
        "step",
        "stepsrc",
    }

    @property
    def color(self):
        """
        Sets the color for each cluster step.

        The 'color' property is a color and may be specified as:
          - A hex string (e.g. '#ff0000')
          - An rgb/rgba string (e.g. 'rgb(255,0,0)')
          - An hsl/hsla string (e.g. 'hsl(0,100%,50%)')
          - An hsv/hsva string (e.g. 'hsv(0,100%,100%)')
          - A named CSS color: see https://plotly.com/python/css-colors/ for a list
          - A list or array of any of the above

        Returns
        -------
        str|numpy.ndarray
        """
        return self["color"]

    @color.setter
    def color(self, val):
        self["color"] = val

    @property
    def colorsrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `color`.

        The 'colorsrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["colorsrc"]

    @colorsrc.setter
    def colorsrc(self, val):
        self["colorsrc"] = val

    @property
    def enabled(self):
        """
        Determines whether clustering is enabled or disabled.

        The 'enabled' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["enabled"]

    @enabled.setter
    def enabled(self, val):
        self["enabled"] = val

    @property
    def maxzoom(self):
        """
        Sets the maximum zoom level. At zoom levels equal to or greater
        than this, points will never be clustered.

        The 'maxzoom' property is a number and may be specified as:
          - An int or float in the interval [0, 24]

        Returns
        -------
        int|float
        """
        return self["maxzoom"]

    @maxzoom.setter
    def maxzoom(self, val):
        self["maxzoom"] = val

    @property
    def opacity(self):
        """
        Sets the marker opacity.

        The 'opacity' property is a number and may be specified as:
          - An int or float in the interval [0, 1]
          - A tuple, list, or one-dimensional numpy array of the above

        Returns
        -------
        int|float|numpy.ndarray
        """
        return self["opacity"]

    @opacity.setter
    def opacity(self, val):
        self["opacity"] = val

    @property
    def opacitysrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `opacity`.

        The 'opacitysrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["opacitysrc"]

    @opacitysrc.setter
    def opacitysrc(self, val):
        self["opacitysrc"] = val

    @property
    def size(self):
        """
        Sets the size for each cluster step.

        The 'size' property is a number and may be specified as:
          - An int or float in the interval [0, inf]
          - A tuple, list, or one-dimensional numpy array of the above

        Returns
        -------
        int|float|numpy.ndarray
        """
        return self["size"]

    @size.setter
    def size(self, val):
        self["size"] = val

    @property
    def sizesrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `size`.

        The 'sizesrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["sizesrc"]

    @sizesrc.setter
    def sizesrc(self, val):
        self["sizesrc"] = val

    @property
    def step(self):
        """
        Sets how many points it takes to create a cluster or advance to
        the next cluster step. Use this in conjunction with arrays for
        `size` and / or `color`. If an integer, steps start at
        multiples of this number. If an array, each step extends from
        the given value until one less than the next value.

        The 'step' property is a number and may be specified as:
          - An int or float in the interval [-1, inf]
          - A tuple, list, or one-dimensional numpy array of the above

        Returns
        -------
        int|float|numpy.ndarray
        """
        return self["step"]

    @step.setter
    def step(self, val):
        self["step"] = val

    @property
    def stepsrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `step`.

        The 'stepsrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["stepsrc"]

    @stepsrc.setter
    def stepsrc(self, val):
        self["stepsrc"] = val

    @property
    def _prop_descriptions(self):
        return """\
        color
            Sets the color for each cluster step.
        colorsrc
            Sets the source reference on Chart Studio Cloud for
            `color`.
        enabled
            Determines whether clustering is enabled or disabled.
        maxzoom
            Sets the maximum zoom level. At zoom levels equal to or
            greater than this, points will never be clustered.
        opacity
            Sets the marker opacity.
        opacitysrc
            Sets the source reference on Chart Studio Cloud for
            `opacity`.
        size
            Sets the size for each cluster step.
        sizesrc
            Sets the source reference on Chart Studio Cloud for
            `size`.
        step
            Sets how many points it takes to create a cluster or
            advance to the next cluster step. Use this in
            conjunction with arrays for `size` and / or `color`. If
            an integer, steps start at multiples of this number. If
            an array, each step extends from the given value until
            one less than the next value.
        stepsrc
            Sets the source reference on Chart Studio Cloud for
            `step`.
        """

    def __init__(
        self,
        arg=None,
        color=None,
        colorsrc=None,
        enabled=None,
        maxzoom=None,
        opacity=None,
        opacitysrc=None,
        size=None,
        sizesrc=None,
        step=None,
        stepsrc=None,
        **kwargs,
    ):
        """
        Construct a new Cluster object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.scattermap.Cluster`
        color
            Sets the color for each cluster step.
        colorsrc
            Sets the source reference on Chart Studio Cloud for
            `color`.
        enabled
            Determines whether clustering is enabled or disabled.
        maxzoom
            Sets the maximum zoom level. At zoom levels equal to or
            greater than this, points will never be clustered.
        opacity
            Sets the marker opacity.
        opacitysrc
            Sets the source reference on Chart Studio Cloud for
            `opacity`.
        size
            Sets the size for each cluster step.
        sizesrc
            Sets the source reference on Chart Studio Cloud for
            `size`.
        step
            Sets how many points it takes to create a cluster or
            advance to the next cluster step. Use this in
            conjunction with arrays for `size` and / or `color`. If
            an integer, steps start at multiples of this number. If
            an array, each step extends from the given value until
            one less than the next value.
        stepsrc
            Sets the source reference on Chart Studio Cloud for
            `step`.

        Returns
        -------
        Cluster
        """
        super().__init__("cluster")
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
The first argument to the plotly.graph_objs.scattermap.Cluster
constructor must be a dict or
an instance of :class:`plotly.graph_objs.scattermap.Cluster`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("color", arg, color)
        self._set_property("colorsrc", arg, colorsrc)
        self._set_property("enabled", arg, enabled)
        self._set_property("maxzoom", arg, maxzoom)
        self._set_property("opacity", arg, opacity)
        self._set_property("opacitysrc", arg, opacitysrc)
        self._set_property("size", arg, size)
        self._set_property("sizesrc", arg, sizesrc)
        self._set_property("step", arg, step)
        self._set_property("stepsrc", arg, stepsrc)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
