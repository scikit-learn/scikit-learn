from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Marker(_BaseTraceHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "pointcloud"
    _path_str = "pointcloud.marker"
    _valid_props = {"blend", "border", "color", "opacity", "sizemax", "sizemin"}

    # blend
    # -----
    @property
    def blend(self):
        """
        Determines if colors are blended together for a translucency
        effect in case `opacity` is specified as a value less then `1`.
        Setting `blend` to `true` reduces zoom/pan speed if used with
        large numbers of points.

        The 'blend' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["blend"]

    @blend.setter
    def blend(self, val):
        self["blend"] = val

    # border
    # ------
    @property
    def border(self):
        """
        The 'border' property is an instance of Border
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.pointcloud.marker.Border`
          - A dict of string/value properties that will be passed
            to the Border constructor

            Supported dict properties:

                arearatio
                    Specifies what fraction of the marker area is
                    covered with the border.
                color
                    Sets the stroke color. It accepts a specific
                    color. If the color is not fully opaque and
                    there are hundreds of thousands of points, it
                    may cause slower zooming and panning.

        Returns
        -------
        plotly.graph_objs.pointcloud.marker.Border
        """
        return self["border"]

    @border.setter
    def border(self, val):
        self["border"] = val

    # color
    # -----
    @property
    def color(self):
        """
        Sets the marker fill color. It accepts a specific color. If the
        color is not fully opaque and there are hundreds of thousands
        of points, it may cause slower zooming and panning.

        The 'color' property is a color and may be specified as:
          - A hex string (e.g. '#ff0000')
          - An rgb/rgba string (e.g. 'rgb(255,0,0)')
          - An hsl/hsla string (e.g. 'hsl(0,100%,50%)')
          - An hsv/hsva string (e.g. 'hsv(0,100%,100%)')
          - A named CSS color:
                aliceblue, antiquewhite, aqua, aquamarine, azure,
                beige, bisque, black, blanchedalmond, blue,
                blueviolet, brown, burlywood, cadetblue,
                chartreuse, chocolate, coral, cornflowerblue,
                cornsilk, crimson, cyan, darkblue, darkcyan,
                darkgoldenrod, darkgray, darkgrey, darkgreen,
                darkkhaki, darkmagenta, darkolivegreen, darkorange,
                darkorchid, darkred, darksalmon, darkseagreen,
                darkslateblue, darkslategray, darkslategrey,
                darkturquoise, darkviolet, deeppink, deepskyblue,
                dimgray, dimgrey, dodgerblue, firebrick,
                floralwhite, forestgreen, fuchsia, gainsboro,
                ghostwhite, gold, goldenrod, gray, grey, green,
                greenyellow, honeydew, hotpink, indianred, indigo,
                ivory, khaki, lavender, lavenderblush, lawngreen,
                lemonchiffon, lightblue, lightcoral, lightcyan,
                lightgoldenrodyellow, lightgray, lightgrey,
                lightgreen, lightpink, lightsalmon, lightseagreen,
                lightskyblue, lightslategray, lightslategrey,
                lightsteelblue, lightyellow, lime, limegreen,
                linen, magenta, maroon, mediumaquamarine,
                mediumblue, mediumorchid, mediumpurple,
                mediumseagreen, mediumslateblue, mediumspringgreen,
                mediumturquoise, mediumvioletred, midnightblue,
                mintcream, mistyrose, moccasin, navajowhite, navy,
                oldlace, olive, olivedrab, orange, orangered,
                orchid, palegoldenrod, palegreen, paleturquoise,
                palevioletred, papayawhip, peachpuff, peru, pink,
                plum, powderblue, purple, red, rosybrown,
                royalblue, rebeccapurple, saddlebrown, salmon,
                sandybrown, seagreen, seashell, sienna, silver,
                skyblue, slateblue, slategray, slategrey, snow,
                springgreen, steelblue, tan, teal, thistle, tomato,
                turquoise, violet, wheat, white, whitesmoke,
                yellow, yellowgreen

        Returns
        -------
        str
        """
        return self["color"]

    @color.setter
    def color(self, val):
        self["color"] = val

    # opacity
    # -------
    @property
    def opacity(self):
        """
        Sets the marker opacity. The default value is `1` (fully
        opaque). If the markers are not fully opaque and there are
        hundreds of thousands of points, it may cause slower zooming
        and panning. Opacity fades the color even if `blend` is left on
        `false` even if there is no translucency effect in that case.

        The 'opacity' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["opacity"]

    @opacity.setter
    def opacity(self, val):
        self["opacity"] = val

    # sizemax
    # -------
    @property
    def sizemax(self):
        """
        Sets the maximum size (in px) of the rendered marker points.
        Effective when the `pointcloud` shows only few points.

        The 'sizemax' property is a number and may be specified as:
          - An int or float in the interval [0.1, inf]

        Returns
        -------
        int|float
        """
        return self["sizemax"]

    @sizemax.setter
    def sizemax(self, val):
        self["sizemax"] = val

    # sizemin
    # -------
    @property
    def sizemin(self):
        """
        Sets the minimum size (in px) of the rendered marker points,
        effective when the `pointcloud` shows a million or more points.

        The 'sizemin' property is a number and may be specified as:
          - An int or float in the interval [0.1, 2]

        Returns
        -------
        int|float
        """
        return self["sizemin"]

    @sizemin.setter
    def sizemin(self, val):
        self["sizemin"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        blend
            Determines if colors are blended together for a
            translucency effect in case `opacity` is specified as a
            value less then `1`. Setting `blend` to `true` reduces
            zoom/pan speed if used with large numbers of points.
        border
            :class:`plotly.graph_objects.pointcloud.marker.Border`
            instance or dict with compatible properties
        color
            Sets the marker fill color. It accepts a specific
            color. If the color is not fully opaque and there are
            hundreds of thousands of points, it may cause slower
            zooming and panning.
        opacity
            Sets the marker opacity. The default value is `1`
            (fully opaque). If the markers are not fully opaque and
            there are hundreds of thousands of points, it may cause
            slower zooming and panning. Opacity fades the color
            even if `blend` is left on `false` even if there is no
            translucency effect in that case.
        sizemax
            Sets the maximum size (in px) of the rendered marker
            points. Effective when the `pointcloud` shows only few
            points.
        sizemin
            Sets the minimum size (in px) of the rendered marker
            points, effective when the `pointcloud` shows a million
            or more points.
        """

    def __init__(
        self,
        arg=None,
        blend=None,
        border=None,
        color=None,
        opacity=None,
        sizemax=None,
        sizemin=None,
        **kwargs,
    ):
        """
        Construct a new Marker object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.pointcloud.Marker`
        blend
            Determines if colors are blended together for a
            translucency effect in case `opacity` is specified as a
            value less then `1`. Setting `blend` to `true` reduces
            zoom/pan speed if used with large numbers of points.
        border
            :class:`plotly.graph_objects.pointcloud.marker.Border`
            instance or dict with compatible properties
        color
            Sets the marker fill color. It accepts a specific
            color. If the color is not fully opaque and there are
            hundreds of thousands of points, it may cause slower
            zooming and panning.
        opacity
            Sets the marker opacity. The default value is `1`
            (fully opaque). If the markers are not fully opaque and
            there are hundreds of thousands of points, it may cause
            slower zooming and panning. Opacity fades the color
            even if `blend` is left on `false` even if there is no
            translucency effect in that case.
        sizemax
            Sets the maximum size (in px) of the rendered marker
            points. Effective when the `pointcloud` shows only few
            points.
        sizemin
            Sets the minimum size (in px) of the rendered marker
            points, effective when the `pointcloud` shows a million
            or more points.

        Returns
        -------
        Marker
        """
        super(Marker, self).__init__("marker")

        if "_parent" in kwargs:
            self._parent = kwargs["_parent"]
            return

        # Validate arg
        # ------------
        if arg is None:
            arg = {}
        elif isinstance(arg, self.__class__):
            arg = arg.to_plotly_json()
        elif isinstance(arg, dict):
            arg = _copy.copy(arg)
        else:
            raise ValueError(
                """\
The first argument to the plotly.graph_objs.pointcloud.Marker
constructor must be a dict or
an instance of :class:`plotly.graph_objs.pointcloud.Marker`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("blend", None)
        _v = blend if blend is not None else _v
        if _v is not None:
            self["blend"] = _v
        _v = arg.pop("border", None)
        _v = border if border is not None else _v
        if _v is not None:
            self["border"] = _v
        _v = arg.pop("color", None)
        _v = color if color is not None else _v
        if _v is not None:
            self["color"] = _v
        _v = arg.pop("opacity", None)
        _v = opacity if opacity is not None else _v
        if _v is not None:
            self["opacity"] = _v
        _v = arg.pop("sizemax", None)
        _v = sizemax if sizemax is not None else _v
        if _v is not None:
            self["sizemax"] = _v
        _v = arg.pop("sizemin", None)
        _v = sizemin if sizemin is not None else _v
        if _v is not None:
            self["sizemin"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
