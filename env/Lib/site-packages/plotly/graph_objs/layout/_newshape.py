from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Newshape(_BaseLayoutHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "layout"
    _path_str = "layout.newshape"
    _valid_props = {
        "drawdirection",
        "fillcolor",
        "fillrule",
        "label",
        "layer",
        "legend",
        "legendgroup",
        "legendgrouptitle",
        "legendrank",
        "legendwidth",
        "line",
        "name",
        "opacity",
        "showlegend",
        "visible",
    }

    # drawdirection
    # -------------
    @property
    def drawdirection(self):
        """
        When `dragmode` is set to "drawrect", "drawline" or
        "drawcircle" this limits the drag to be horizontal, vertical or
        diagonal. Using "diagonal" there is no limit e.g. in drawing
        lines in any direction. "ortho" limits the draw to be either
        horizontal or vertical. "horizontal" allows horizontal extend.
        "vertical" allows vertical extend.

        The 'drawdirection' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['ortho', 'horizontal', 'vertical', 'diagonal']

        Returns
        -------
        Any
        """
        return self["drawdirection"]

    @drawdirection.setter
    def drawdirection(self, val):
        self["drawdirection"] = val

    # fillcolor
    # ---------
    @property
    def fillcolor(self):
        """
        Sets the color filling new shapes' interior. Please note that
        if using a fillcolor with alpha greater than half, drag inside
        the active shape starts moving the shape underneath, otherwise
        a new shape could be started over.

        The 'fillcolor' property is a color and may be specified as:
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
        return self["fillcolor"]

    @fillcolor.setter
    def fillcolor(self, val):
        self["fillcolor"] = val

    # fillrule
    # --------
    @property
    def fillrule(self):
        """
        Determines the path's interior. For more info please visit
        https://developer.mozilla.org/en-
        US/docs/Web/SVG/Attribute/fill-rule

        The 'fillrule' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['evenodd', 'nonzero']

        Returns
        -------
        Any
        """
        return self["fillrule"]

    @fillrule.setter
    def fillrule(self, val):
        self["fillrule"] = val

    # label
    # -----
    @property
    def label(self):
        """
        The 'label' property is an instance of Label
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.newshape.Label`
          - A dict of string/value properties that will be passed
            to the Label constructor

            Supported dict properties:

                font
                    Sets the new shape label text font.
                padding
                    Sets padding (in px) between edge of label and
                    edge of new shape.
                text
                    Sets the text to display with the new shape. It
                    is also used for legend item if `name` is not
                    provided.
                textangle
                    Sets the angle at which the label text is drawn
                    with respect to the horizontal. For lines,
                    angle "auto" is the same angle as the line. For
                    all other shapes, angle "auto" is horizontal.
                textposition
                    Sets the position of the label text relative to
                    the new shape. Supported values for rectangles,
                    circles and paths are *top left*, *top center*,
                    *top right*, *middle left*, *middle center*,
                    *middle right*, *bottom left*, *bottom center*,
                    and *bottom right*. Supported values for lines
                    are "start", "middle", and "end". Default:
                    *middle center* for rectangles, circles, and
                    paths; "middle" for lines.
                texttemplate
                    Template string used for rendering the new
                    shape's label. Note that this will override
                    `text`. Variables are inserted using
                    %{variable}, for example "x0: %{x0}". Numbers
                    are formatted using d3-format's syntax
                    %{variable:d3-format}, for example "Price:
                    %{x0:$.2f}". See https://github.com/d3/d3-
                    format/tree/v1.4.5#d3-format for details on the
                    formatting syntax. Dates are formatted using
                    d3-time-format's syntax %{variable|d3-time-
                    format}, for example "Day: %{x0|%m %b %Y}". See
                    https://github.com/d3/d3-time-
                    format/tree/v2.2.3#locale_format for details on
                    the date formatting syntax. A single
                    multiplication or division operation may be
                    applied to numeric variables, and combined with
                    d3 number formatting, for example "Length in
                    cm: %{x0*2.54}", "%{slope*60:.1f} meters per
                    second." For log axes, variable values are
                    given in log units. For date axes, x/y
                    coordinate variables and center variables use
                    datetimes, while all other variable values use
                    values in ms. Finally, the template string has
                    access to variables `x0`, `x1`, `y0`, `y1`,
                    `slope`, `dx`, `dy`, `width`, `height`,
                    `length`, `xcenter` and `ycenter`.
                xanchor
                    Sets the label's horizontal position anchor
                    This anchor binds the specified `textposition`
                    to the "left", "center" or "right" of the label
                    text. For example, if `textposition` is set to
                    *top right* and `xanchor` to "right" then the
                    right-most portion of the label text lines up
                    with the right-most edge of the new shape.
                yanchor
                    Sets the label's vertical position anchor This
                    anchor binds the specified `textposition` to
                    the "top", "middle" or "bottom" of the label
                    text. For example, if `textposition` is set to
                    *top right* and `yanchor` to "top" then the
                    top-most portion of the label text lines up
                    with the top-most edge of the new shape.

        Returns
        -------
        plotly.graph_objs.layout.newshape.Label
        """
        return self["label"]

    @label.setter
    def label(self, val):
        self["label"] = val

    # layer
    # -----
    @property
    def layer(self):
        """
        Specifies whether new shapes are drawn below gridlines
        ("below"), between gridlines and traces ("between") or above
        traces ("above").

        The 'layer' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['below', 'above', 'between']

        Returns
        -------
        Any
        """
        return self["layer"]

    @layer.setter
    def layer(self, val):
        self["layer"] = val

    # legend
    # ------
    @property
    def legend(self):
        """
        Sets the reference to a legend to show new shape in. References
        to these legends are "legend", "legend2", "legend3", etc.
        Settings for these legends are set in the layout, under
        `layout.legend`, `layout.legend2`, etc.

        The 'legend' property is an identifier of a particular
        subplot, of type 'legend', that may be specified as the string 'legend'
        optionally followed by an integer >= 1
        (e.g. 'legend', 'legend1', 'legend2', 'legend3', etc.)

        Returns
        -------
        str
        """
        return self["legend"]

    @legend.setter
    def legend(self, val):
        self["legend"] = val

    # legendgroup
    # -----------
    @property
    def legendgroup(self):
        """
        Sets the legend group for new shape. Traces and shapes part of
        the same legend group hide/show at the same time when toggling
        legend items.

        The 'legendgroup' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["legendgroup"]

    @legendgroup.setter
    def legendgroup(self, val):
        self["legendgroup"] = val

    # legendgrouptitle
    # ----------------
    @property
    def legendgrouptitle(self):
        """
        The 'legendgrouptitle' property is an instance of Legendgrouptitle
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.newshape.Legendgrouptitle`
          - A dict of string/value properties that will be passed
            to the Legendgrouptitle constructor

            Supported dict properties:

                font
                    Sets this legend group's title font.
                text
                    Sets the title of the legend group.

        Returns
        -------
        plotly.graph_objs.layout.newshape.Legendgrouptitle
        """
        return self["legendgrouptitle"]

    @legendgrouptitle.setter
    def legendgrouptitle(self, val):
        self["legendgrouptitle"] = val

    # legendrank
    # ----------
    @property
    def legendrank(self):
        """
        Sets the legend rank for new shape. Items and groups with
        smaller ranks are presented on top/left side while with
        "reversed" `legend.traceorder` they are on bottom/right side.
        The default legendrank is 1000, so that you can use ranks less
        than 1000 to place certain items before all unranked items, and
        ranks greater than 1000 to go after all unranked items.

        The 'legendrank' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["legendrank"]

    @legendrank.setter
    def legendrank(self, val):
        self["legendrank"] = val

    # legendwidth
    # -----------
    @property
    def legendwidth(self):
        """
        Sets the width (in px or fraction) of the legend for new shape.

        The 'legendwidth' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["legendwidth"]

    @legendwidth.setter
    def legendwidth(self, val):
        self["legendwidth"] = val

    # line
    # ----
    @property
    def line(self):
        """
        The 'line' property is an instance of Line
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.newshape.Line`
          - A dict of string/value properties that will be passed
            to the Line constructor

            Supported dict properties:

                color
                    Sets the line color. By default uses either
                    dark grey or white to increase contrast with
                    background color.
                dash
                    Sets the dash style of lines. Set to a dash
                    type string ("solid", "dot", "dash",
                    "longdash", "dashdot", or "longdashdot") or a
                    dash length list in px (eg "5px,10px,2px,2px").
                width
                    Sets the line width (in px).

        Returns
        -------
        plotly.graph_objs.layout.newshape.Line
        """
        return self["line"]

    @line.setter
    def line(self, val):
        self["line"] = val

    # name
    # ----
    @property
    def name(self):
        """
        Sets new shape name. The name appears as the legend item.

        The 'name' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["name"]

    @name.setter
    def name(self, val):
        self["name"] = val

    # opacity
    # -------
    @property
    def opacity(self):
        """
        Sets the opacity of new shapes.

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

    # showlegend
    # ----------
    @property
    def showlegend(self):
        """
        Determines whether or not new shape is shown in the legend.

        The 'showlegend' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["showlegend"]

    @showlegend.setter
    def showlegend(self, val):
        self["showlegend"] = val

    # visible
    # -------
    @property
    def visible(self):
        """
        Determines whether or not new shape is visible. If
        "legendonly", the shape is not drawn, but can appear as a
        legend item (provided that the legend itself is visible).

        The 'visible' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                [True, False, 'legendonly']

        Returns
        -------
        Any
        """
        return self["visible"]

    @visible.setter
    def visible(self, val):
        self["visible"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        drawdirection
            When `dragmode` is set to "drawrect", "drawline" or
            "drawcircle" this limits the drag to be horizontal,
            vertical or diagonal. Using "diagonal" there is no
            limit e.g. in drawing lines in any direction. "ortho"
            limits the draw to be either horizontal or vertical.
            "horizontal" allows horizontal extend. "vertical"
            allows vertical extend.
        fillcolor
            Sets the color filling new shapes' interior. Please
            note that if using a fillcolor with alpha greater than
            half, drag inside the active shape starts moving the
            shape underneath, otherwise a new shape could be
            started over.
        fillrule
            Determines the path's interior. For more info please
            visit https://developer.mozilla.org/en-
            US/docs/Web/SVG/Attribute/fill-rule
        label
            :class:`plotly.graph_objects.layout.newshape.Label`
            instance or dict with compatible properties
        layer
            Specifies whether new shapes are drawn below gridlines
            ("below"), between gridlines and traces ("between") or
            above traces ("above").
        legend
            Sets the reference to a legend to show new shape in.
            References to these legends are "legend", "legend2",
            "legend3", etc. Settings for these legends are set in
            the layout, under `layout.legend`, `layout.legend2`,
            etc.
        legendgroup
            Sets the legend group for new shape. Traces and shapes
            part of the same legend group hide/show at the same
            time when toggling legend items.
        legendgrouptitle
            :class:`plotly.graph_objects.layout.newshape.Legendgrou
            ptitle` instance or dict with compatible properties
        legendrank
            Sets the legend rank for new shape. Items and groups
            with smaller ranks are presented on top/left side while
            with "reversed" `legend.traceorder` they are on
            bottom/right side. The default legendrank is 1000, so
            that you can use ranks less than 1000 to place certain
            items before all unranked items, and ranks greater than
            1000 to go after all unranked items.
        legendwidth
            Sets the width (in px or fraction) of the legend for
            new shape.
        line
            :class:`plotly.graph_objects.layout.newshape.Line`
            instance or dict with compatible properties
        name
            Sets new shape name. The name appears as the legend
            item.
        opacity
            Sets the opacity of new shapes.
        showlegend
            Determines whether or not new shape is shown in the
            legend.
        visible
            Determines whether or not new shape is visible. If
            "legendonly", the shape is not drawn, but can appear as
            a legend item (provided that the legend itself is
            visible).
        """

    def __init__(
        self,
        arg=None,
        drawdirection=None,
        fillcolor=None,
        fillrule=None,
        label=None,
        layer=None,
        legend=None,
        legendgroup=None,
        legendgrouptitle=None,
        legendrank=None,
        legendwidth=None,
        line=None,
        name=None,
        opacity=None,
        showlegend=None,
        visible=None,
        **kwargs,
    ):
        """
        Construct a new Newshape object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.layout.Newshape`
        drawdirection
            When `dragmode` is set to "drawrect", "drawline" or
            "drawcircle" this limits the drag to be horizontal,
            vertical or diagonal. Using "diagonal" there is no
            limit e.g. in drawing lines in any direction. "ortho"
            limits the draw to be either horizontal or vertical.
            "horizontal" allows horizontal extend. "vertical"
            allows vertical extend.
        fillcolor
            Sets the color filling new shapes' interior. Please
            note that if using a fillcolor with alpha greater than
            half, drag inside the active shape starts moving the
            shape underneath, otherwise a new shape could be
            started over.
        fillrule
            Determines the path's interior. For more info please
            visit https://developer.mozilla.org/en-
            US/docs/Web/SVG/Attribute/fill-rule
        label
            :class:`plotly.graph_objects.layout.newshape.Label`
            instance or dict with compatible properties
        layer
            Specifies whether new shapes are drawn below gridlines
            ("below"), between gridlines and traces ("between") or
            above traces ("above").
        legend
            Sets the reference to a legend to show new shape in.
            References to these legends are "legend", "legend2",
            "legend3", etc. Settings for these legends are set in
            the layout, under `layout.legend`, `layout.legend2`,
            etc.
        legendgroup
            Sets the legend group for new shape. Traces and shapes
            part of the same legend group hide/show at the same
            time when toggling legend items.
        legendgrouptitle
            :class:`plotly.graph_objects.layout.newshape.Legendgrou
            ptitle` instance or dict with compatible properties
        legendrank
            Sets the legend rank for new shape. Items and groups
            with smaller ranks are presented on top/left side while
            with "reversed" `legend.traceorder` they are on
            bottom/right side. The default legendrank is 1000, so
            that you can use ranks less than 1000 to place certain
            items before all unranked items, and ranks greater than
            1000 to go after all unranked items.
        legendwidth
            Sets the width (in px or fraction) of the legend for
            new shape.
        line
            :class:`plotly.graph_objects.layout.newshape.Line`
            instance or dict with compatible properties
        name
            Sets new shape name. The name appears as the legend
            item.
        opacity
            Sets the opacity of new shapes.
        showlegend
            Determines whether or not new shape is shown in the
            legend.
        visible
            Determines whether or not new shape is visible. If
            "legendonly", the shape is not drawn, but can appear as
            a legend item (provided that the legend itself is
            visible).

        Returns
        -------
        Newshape
        """
        super(Newshape, self).__init__("newshape")

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
The first argument to the plotly.graph_objs.layout.Newshape
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.Newshape`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("drawdirection", None)
        _v = drawdirection if drawdirection is not None else _v
        if _v is not None:
            self["drawdirection"] = _v
        _v = arg.pop("fillcolor", None)
        _v = fillcolor if fillcolor is not None else _v
        if _v is not None:
            self["fillcolor"] = _v
        _v = arg.pop("fillrule", None)
        _v = fillrule if fillrule is not None else _v
        if _v is not None:
            self["fillrule"] = _v
        _v = arg.pop("label", None)
        _v = label if label is not None else _v
        if _v is not None:
            self["label"] = _v
        _v = arg.pop("layer", None)
        _v = layer if layer is not None else _v
        if _v is not None:
            self["layer"] = _v
        _v = arg.pop("legend", None)
        _v = legend if legend is not None else _v
        if _v is not None:
            self["legend"] = _v
        _v = arg.pop("legendgroup", None)
        _v = legendgroup if legendgroup is not None else _v
        if _v is not None:
            self["legendgroup"] = _v
        _v = arg.pop("legendgrouptitle", None)
        _v = legendgrouptitle if legendgrouptitle is not None else _v
        if _v is not None:
            self["legendgrouptitle"] = _v
        _v = arg.pop("legendrank", None)
        _v = legendrank if legendrank is not None else _v
        if _v is not None:
            self["legendrank"] = _v
        _v = arg.pop("legendwidth", None)
        _v = legendwidth if legendwidth is not None else _v
        if _v is not None:
            self["legendwidth"] = _v
        _v = arg.pop("line", None)
        _v = line if line is not None else _v
        if _v is not None:
            self["line"] = _v
        _v = arg.pop("name", None)
        _v = name if name is not None else _v
        if _v is not None:
            self["name"] = _v
        _v = arg.pop("opacity", None)
        _v = opacity if opacity is not None else _v
        if _v is not None:
            self["opacity"] = _v
        _v = arg.pop("showlegend", None)
        _v = showlegend if showlegend is not None else _v
        if _v is not None:
            self["showlegend"] = _v
        _v = arg.pop("visible", None)
        _v = visible if visible is not None else _v
        if _v is not None:
            self["visible"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
