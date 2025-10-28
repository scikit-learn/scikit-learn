#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Modebar(_BaseLayoutHierarchyType):
    _parent_path_str = "layout"
    _path_str = "layout.modebar"
    _valid_props = {
        "activecolor",
        "add",
        "addsrc",
        "bgcolor",
        "color",
        "orientation",
        "remove",
        "removesrc",
        "uirevision",
    }

    @property
    def activecolor(self):
        """
        Sets the color of the active or hovered on icons in the
        modebar.

        The 'activecolor' property is a color and may be specified as:
          - A hex string (e.g. '#ff0000')
          - An rgb/rgba string (e.g. 'rgb(255,0,0)')
          - An hsl/hsla string (e.g. 'hsl(0,100%,50%)')
          - An hsv/hsva string (e.g. 'hsv(0,100%,100%)')
          - A named CSS color: see https://plotly.com/python/css-colors/ for a list

        Returns
        -------
        str
        """
        return self["activecolor"]

    @activecolor.setter
    def activecolor(self, val):
        self["activecolor"] = val

    @property
    def add(self):
        """
        Determines which predefined modebar buttons to add. Please note
        that these buttons will only be shown if they are compatible
        with all trace types used in a graph. Similar to
        `config.modeBarButtonsToAdd` option. This may include
        "v1hovermode", "hoverclosest", "hovercompare", "togglehover",
        "togglespikelines", "drawline", "drawopenpath",
        "drawclosedpath", "drawcircle", "drawrect", "eraseshape".

        The 'add' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string
          - A tuple, list, or one-dimensional numpy array of the above

        Returns
        -------
        str|numpy.ndarray
        """
        return self["add"]

    @add.setter
    def add(self, val):
        self["add"] = val

    @property
    def addsrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `add`.

        The 'addsrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["addsrc"]

    @addsrc.setter
    def addsrc(self, val):
        self["addsrc"] = val

    @property
    def bgcolor(self):
        """
        Sets the background color of the modebar.

        The 'bgcolor' property is a color and may be specified as:
          - A hex string (e.g. '#ff0000')
          - An rgb/rgba string (e.g. 'rgb(255,0,0)')
          - An hsl/hsla string (e.g. 'hsl(0,100%,50%)')
          - An hsv/hsva string (e.g. 'hsv(0,100%,100%)')
          - A named CSS color: see https://plotly.com/python/css-colors/ for a list

        Returns
        -------
        str
        """
        return self["bgcolor"]

    @bgcolor.setter
    def bgcolor(self, val):
        self["bgcolor"] = val

    @property
    def color(self):
        """
        Sets the color of the icons in the modebar.

        The 'color' property is a color and may be specified as:
          - A hex string (e.g. '#ff0000')
          - An rgb/rgba string (e.g. 'rgb(255,0,0)')
          - An hsl/hsla string (e.g. 'hsl(0,100%,50%)')
          - An hsv/hsva string (e.g. 'hsv(0,100%,100%)')
          - A named CSS color: see https://plotly.com/python/css-colors/ for a list

        Returns
        -------
        str
        """
        return self["color"]

    @color.setter
    def color(self, val):
        self["color"] = val

    @property
    def orientation(self):
        """
        Sets the orientation of the modebar.

        The 'orientation' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['v', 'h']

        Returns
        -------
        Any
        """
        return self["orientation"]

    @orientation.setter
    def orientation(self, val):
        self["orientation"] = val

    @property
    def remove(self):
        """
        Determines which predefined modebar buttons to remove. Similar
        to `config.modeBarButtonsToRemove` option. This may include
        "autoScale2d", "autoscale", "editInChartStudio",
        "editinchartstudio", "hoverCompareCartesian", "hovercompare",
        "lasso", "lasso2d", "orbitRotation", "orbitrotation", "pan",
        "pan2d", "pan3d", "reset", "resetCameraDefault3d",
        "resetCameraLastSave3d", "resetGeo", "resetSankeyGroup",
        "resetScale2d", "resetViewMap", "resetViewMapbox",
        "resetViews", "resetcameradefault", "resetcameralastsave",
        "resetsankeygroup", "resetscale", "resetview", "resetviews",
        "select", "select2d", "sendDataToCloud", "senddatatocloud",
        "tableRotation", "tablerotation", "toImage", "toggleHover",
        "toggleSpikelines", "togglehover", "togglespikelines",
        "toimage", "zoom", "zoom2d", "zoom3d", "zoomIn2d", "zoomInGeo",
        "zoomInMap", "zoomInMapbox", "zoomOut2d", "zoomOutGeo",
        "zoomOutMap", "zoomOutMapbox", "zoomin", "zoomout".

        The 'remove' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string
          - A tuple, list, or one-dimensional numpy array of the above

        Returns
        -------
        str|numpy.ndarray
        """
        return self["remove"]

    @remove.setter
    def remove(self, val):
        self["remove"] = val

    @property
    def removesrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `remove`.

        The 'removesrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["removesrc"]

    @removesrc.setter
    def removesrc(self, val):
        self["removesrc"] = val

    @property
    def uirevision(self):
        """
        Controls persistence of user-driven changes related to the
        modebar, including `hovermode`, `dragmode`, and `showspikes` at
        both the root level and inside subplots. Defaults to
        `layout.uirevision`.

        The 'uirevision' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["uirevision"]

    @uirevision.setter
    def uirevision(self, val):
        self["uirevision"] = val

    @property
    def _prop_descriptions(self):
        return """\
        activecolor
            Sets the color of the active or hovered on icons in the
            modebar.
        add
            Determines which predefined modebar buttons to add.
            Please note that these buttons will only be shown if
            they are compatible with all trace types used in a
            graph. Similar to `config.modeBarButtonsToAdd` option.
            This may include "v1hovermode", "hoverclosest",
            "hovercompare", "togglehover", "togglespikelines",
            "drawline", "drawopenpath", "drawclosedpath",
            "drawcircle", "drawrect", "eraseshape".
        addsrc
            Sets the source reference on Chart Studio Cloud for
            `add`.
        bgcolor
            Sets the background color of the modebar.
        color
            Sets the color of the icons in the modebar.
        orientation
            Sets the orientation of the modebar.
        remove
            Determines which predefined modebar buttons to remove.
            Similar to `config.modeBarButtonsToRemove` option. This
            may include "autoScale2d", "autoscale",
            "editInChartStudio", "editinchartstudio",
            "hoverCompareCartesian", "hovercompare", "lasso",
            "lasso2d", "orbitRotation", "orbitrotation", "pan",
            "pan2d", "pan3d", "reset", "resetCameraDefault3d",
            "resetCameraLastSave3d", "resetGeo",
            "resetSankeyGroup", "resetScale2d", "resetViewMap",
            "resetViewMapbox", "resetViews", "resetcameradefault",
            "resetcameralastsave", "resetsankeygroup",
            "resetscale", "resetview", "resetviews", "select",
            "select2d", "sendDataToCloud", "senddatatocloud",
            "tableRotation", "tablerotation", "toImage",
            "toggleHover", "toggleSpikelines", "togglehover",
            "togglespikelines", "toimage", "zoom", "zoom2d",
            "zoom3d", "zoomIn2d", "zoomInGeo", "zoomInMap",
            "zoomInMapbox", "zoomOut2d", "zoomOutGeo",
            "zoomOutMap", "zoomOutMapbox", "zoomin", "zoomout".
        removesrc
            Sets the source reference on Chart Studio Cloud for
            `remove`.
        uirevision
            Controls persistence of user-driven changes related to
            the modebar, including `hovermode`, `dragmode`, and
            `showspikes` at both the root level and inside
            subplots. Defaults to `layout.uirevision`.
        """

    def __init__(
        self,
        arg=None,
        activecolor=None,
        add=None,
        addsrc=None,
        bgcolor=None,
        color=None,
        orientation=None,
        remove=None,
        removesrc=None,
        uirevision=None,
        **kwargs,
    ):
        """
        Construct a new Modebar object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.layout.Modebar`
        activecolor
            Sets the color of the active or hovered on icons in the
            modebar.
        add
            Determines which predefined modebar buttons to add.
            Please note that these buttons will only be shown if
            they are compatible with all trace types used in a
            graph. Similar to `config.modeBarButtonsToAdd` option.
            This may include "v1hovermode", "hoverclosest",
            "hovercompare", "togglehover", "togglespikelines",
            "drawline", "drawopenpath", "drawclosedpath",
            "drawcircle", "drawrect", "eraseshape".
        addsrc
            Sets the source reference on Chart Studio Cloud for
            `add`.
        bgcolor
            Sets the background color of the modebar.
        color
            Sets the color of the icons in the modebar.
        orientation
            Sets the orientation of the modebar.
        remove
            Determines which predefined modebar buttons to remove.
            Similar to `config.modeBarButtonsToRemove` option. This
            may include "autoScale2d", "autoscale",
            "editInChartStudio", "editinchartstudio",
            "hoverCompareCartesian", "hovercompare", "lasso",
            "lasso2d", "orbitRotation", "orbitrotation", "pan",
            "pan2d", "pan3d", "reset", "resetCameraDefault3d",
            "resetCameraLastSave3d", "resetGeo",
            "resetSankeyGroup", "resetScale2d", "resetViewMap",
            "resetViewMapbox", "resetViews", "resetcameradefault",
            "resetcameralastsave", "resetsankeygroup",
            "resetscale", "resetview", "resetviews", "select",
            "select2d", "sendDataToCloud", "senddatatocloud",
            "tableRotation", "tablerotation", "toImage",
            "toggleHover", "toggleSpikelines", "togglehover",
            "togglespikelines", "toimage", "zoom", "zoom2d",
            "zoom3d", "zoomIn2d", "zoomInGeo", "zoomInMap",
            "zoomInMapbox", "zoomOut2d", "zoomOutGeo",
            "zoomOutMap", "zoomOutMapbox", "zoomin", "zoomout".
        removesrc
            Sets the source reference on Chart Studio Cloud for
            `remove`.
        uirevision
            Controls persistence of user-driven changes related to
            the modebar, including `hovermode`, `dragmode`, and
            `showspikes` at both the root level and inside
            subplots. Defaults to `layout.uirevision`.

        Returns
        -------
        Modebar
        """
        super().__init__("modebar")
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
The first argument to the plotly.graph_objs.layout.Modebar
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.Modebar`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("activecolor", arg, activecolor)
        self._set_property("add", arg, add)
        self._set_property("addsrc", arg, addsrc)
        self._set_property("bgcolor", arg, bgcolor)
        self._set_property("color", arg, color)
        self._set_property("orientation", arg, orientation)
        self._set_property("remove", arg, remove)
        self._set_property("removesrc", arg, removesrc)
        self._set_property("uirevision", arg, uirevision)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
