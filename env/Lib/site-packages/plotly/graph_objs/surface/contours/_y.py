from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Y(_BaseTraceHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "surface.contours"
    _path_str = "surface.contours.y"
    _valid_props = {
        "color",
        "end",
        "highlight",
        "highlightcolor",
        "highlightwidth",
        "project",
        "show",
        "size",
        "start",
        "usecolormap",
        "width",
    }

    # color
    # -----
    @property
    def color(self):
        """
        Sets the color of the contour lines.

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

    # end
    # ---
    @property
    def end(self):
        """
        Sets the end contour level value. Must be more than
        `contours.start`

        The 'end' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["end"]

    @end.setter
    def end(self, val):
        self["end"] = val

    # highlight
    # ---------
    @property
    def highlight(self):
        """
        Determines whether or not contour lines about the y dimension
        are highlighted on hover.

        The 'highlight' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["highlight"]

    @highlight.setter
    def highlight(self, val):
        self["highlight"] = val

    # highlightcolor
    # --------------
    @property
    def highlightcolor(self):
        """
        Sets the color of the highlighted contour lines.

        The 'highlightcolor' property is a color and may be specified as:
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
        return self["highlightcolor"]

    @highlightcolor.setter
    def highlightcolor(self, val):
        self["highlightcolor"] = val

    # highlightwidth
    # --------------
    @property
    def highlightwidth(self):
        """
        Sets the width of the highlighted contour lines.

        The 'highlightwidth' property is a number and may be specified as:
          - An int or float in the interval [1, 16]

        Returns
        -------
        int|float
        """
        return self["highlightwidth"]

    @highlightwidth.setter
    def highlightwidth(self, val):
        self["highlightwidth"] = val

    # project
    # -------
    @property
    def project(self):
        """
        The 'project' property is an instance of Project
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.surface.contours.y.Project`
          - A dict of string/value properties that will be passed
            to the Project constructor

            Supported dict properties:

                x
                    Determines whether or not these contour lines
                    are projected on the x plane. If `highlight` is
                    set to True (the default), the projected lines
                    are shown on hover. If `show` is set to True,
                    the projected lines are shown in permanence.
                y
                    Determines whether or not these contour lines
                    are projected on the y plane. If `highlight` is
                    set to True (the default), the projected lines
                    are shown on hover. If `show` is set to True,
                    the projected lines are shown in permanence.
                z
                    Determines whether or not these contour lines
                    are projected on the z plane. If `highlight` is
                    set to True (the default), the projected lines
                    are shown on hover. If `show` is set to True,
                    the projected lines are shown in permanence.

        Returns
        -------
        plotly.graph_objs.surface.contours.y.Project
        """
        return self["project"]

    @project.setter
    def project(self, val):
        self["project"] = val

    # show
    # ----
    @property
    def show(self):
        """
        Determines whether or not contour lines about the y dimension
        are drawn.

        The 'show' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["show"]

    @show.setter
    def show(self, val):
        self["show"] = val

    # size
    # ----
    @property
    def size(self):
        """
        Sets the step between each contour level. Must be positive.

        The 'size' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["size"]

    @size.setter
    def size(self, val):
        self["size"] = val

    # start
    # -----
    @property
    def start(self):
        """
        Sets the starting contour level value. Must be less than
        `contours.end`

        The 'start' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["start"]

    @start.setter
    def start(self, val):
        self["start"] = val

    # usecolormap
    # -----------
    @property
    def usecolormap(self):
        """
        An alternate to "color". Determines whether or not the contour
        lines are colored using the trace "colorscale".

        The 'usecolormap' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["usecolormap"]

    @usecolormap.setter
    def usecolormap(self, val):
        self["usecolormap"] = val

    # width
    # -----
    @property
    def width(self):
        """
        Sets the width of the contour lines.

        The 'width' property is a number and may be specified as:
          - An int or float in the interval [1, 16]

        Returns
        -------
        int|float
        """
        return self["width"]

    @width.setter
    def width(self, val):
        self["width"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        color
            Sets the color of the contour lines.
        end
            Sets the end contour level value. Must be more than
            `contours.start`
        highlight
            Determines whether or not contour lines about the y
            dimension are highlighted on hover.
        highlightcolor
            Sets the color of the highlighted contour lines.
        highlightwidth
            Sets the width of the highlighted contour lines.
        project
            :class:`plotly.graph_objects.surface.contours.y.Project
            ` instance or dict with compatible properties
        show
            Determines whether or not contour lines about the y
            dimension are drawn.
        size
            Sets the step between each contour level. Must be
            positive.
        start
            Sets the starting contour level value. Must be less
            than `contours.end`
        usecolormap
            An alternate to "color". Determines whether or not the
            contour lines are colored using the trace "colorscale".
        width
            Sets the width of the contour lines.
        """

    def __init__(
        self,
        arg=None,
        color=None,
        end=None,
        highlight=None,
        highlightcolor=None,
        highlightwidth=None,
        project=None,
        show=None,
        size=None,
        start=None,
        usecolormap=None,
        width=None,
        **kwargs,
    ):
        """
        Construct a new Y object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.surface.contours.Y`
        color
            Sets the color of the contour lines.
        end
            Sets the end contour level value. Must be more than
            `contours.start`
        highlight
            Determines whether or not contour lines about the y
            dimension are highlighted on hover.
        highlightcolor
            Sets the color of the highlighted contour lines.
        highlightwidth
            Sets the width of the highlighted contour lines.
        project
            :class:`plotly.graph_objects.surface.contours.y.Project
            ` instance or dict with compatible properties
        show
            Determines whether or not contour lines about the y
            dimension are drawn.
        size
            Sets the step between each contour level. Must be
            positive.
        start
            Sets the starting contour level value. Must be less
            than `contours.end`
        usecolormap
            An alternate to "color". Determines whether or not the
            contour lines are colored using the trace "colorscale".
        width
            Sets the width of the contour lines.

        Returns
        -------
        Y
        """
        super(Y, self).__init__("y")

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
The first argument to the plotly.graph_objs.surface.contours.Y
constructor must be a dict or
an instance of :class:`plotly.graph_objs.surface.contours.Y`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("color", None)
        _v = color if color is not None else _v
        if _v is not None:
            self["color"] = _v
        _v = arg.pop("end", None)
        _v = end if end is not None else _v
        if _v is not None:
            self["end"] = _v
        _v = arg.pop("highlight", None)
        _v = highlight if highlight is not None else _v
        if _v is not None:
            self["highlight"] = _v
        _v = arg.pop("highlightcolor", None)
        _v = highlightcolor if highlightcolor is not None else _v
        if _v is not None:
            self["highlightcolor"] = _v
        _v = arg.pop("highlightwidth", None)
        _v = highlightwidth if highlightwidth is not None else _v
        if _v is not None:
            self["highlightwidth"] = _v
        _v = arg.pop("project", None)
        _v = project if project is not None else _v
        if _v is not None:
            self["project"] = _v
        _v = arg.pop("show", None)
        _v = show if show is not None else _v
        if _v is not None:
            self["show"] = _v
        _v = arg.pop("size", None)
        _v = size if size is not None else _v
        if _v is not None:
            self["size"] = _v
        _v = arg.pop("start", None)
        _v = start if start is not None else _v
        if _v is not None:
            self["start"] = _v
        _v = arg.pop("usecolormap", None)
        _v = usecolormap if usecolormap is not None else _v
        if _v is not None:
            self["usecolormap"] = _v
        _v = arg.pop("width", None)
        _v = width if width is not None else _v
        if _v is not None:
            self["width"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
