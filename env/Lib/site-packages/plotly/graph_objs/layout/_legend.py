from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Legend(_BaseLayoutHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "layout"
    _path_str = "layout.legend"
    _valid_props = {
        "bgcolor",
        "bordercolor",
        "borderwidth",
        "entrywidth",
        "entrywidthmode",
        "font",
        "groupclick",
        "grouptitlefont",
        "indentation",
        "itemclick",
        "itemdoubleclick",
        "itemsizing",
        "itemwidth",
        "orientation",
        "title",
        "tracegroupgap",
        "traceorder",
        "uirevision",
        "valign",
        "visible",
        "x",
        "xanchor",
        "xref",
        "y",
        "yanchor",
        "yref",
    }

    # bgcolor
    # -------
    @property
    def bgcolor(self):
        """
        Sets the legend background color. Defaults to
        `layout.paper_bgcolor`.

        The 'bgcolor' property is a color and may be specified as:
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
        return self["bgcolor"]

    @bgcolor.setter
    def bgcolor(self, val):
        self["bgcolor"] = val

    # bordercolor
    # -----------
    @property
    def bordercolor(self):
        """
        Sets the color of the border enclosing the legend.

        The 'bordercolor' property is a color and may be specified as:
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
        return self["bordercolor"]

    @bordercolor.setter
    def bordercolor(self, val):
        self["bordercolor"] = val

    # borderwidth
    # -----------
    @property
    def borderwidth(self):
        """
        Sets the width (in px) of the border enclosing the legend.

        The 'borderwidth' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["borderwidth"]

    @borderwidth.setter
    def borderwidth(self, val):
        self["borderwidth"] = val

    # entrywidth
    # ----------
    @property
    def entrywidth(self):
        """
        Sets the width (in px or fraction) of the legend. Use 0 to size
        the entry based on the text width, when `entrywidthmode` is set
        to "pixels".

        The 'entrywidth' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["entrywidth"]

    @entrywidth.setter
    def entrywidth(self, val):
        self["entrywidth"] = val

    # entrywidthmode
    # --------------
    @property
    def entrywidthmode(self):
        """
        Determines what entrywidth means.

        The 'entrywidthmode' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['fraction', 'pixels']

        Returns
        -------
        Any
        """
        return self["entrywidthmode"]

    @entrywidthmode.setter
    def entrywidthmode(self, val):
        self["entrywidthmode"] = val

    # font
    # ----
    @property
    def font(self):
        """
        Sets the font used to text the legend items.

        The 'font' property is an instance of Font
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.legend.Font`
          - A dict of string/value properties that will be passed
            to the Font constructor

            Supported dict properties:

                color

                family
                    HTML font family - the typeface that will be
                    applied by the web browser. The web browser
                    will only be able to apply a font if it is
                    available on the system which it operates.
                    Provide multiple font families, separated by
                    commas, to indicate the preference in which to
                    apply fonts if they aren't available on the
                    system. The Chart Studio Cloud (at
                    https://chart-studio.plotly.com or on-premise)
                    generates images on a server, where only a
                    select number of fonts are installed and
                    supported. These include "Arial", "Balto",
                    "Courier New", "Droid Sans", "Droid Serif",
                    "Droid Sans Mono", "Gravitas One", "Old
                    Standard TT", "Open Sans", "Overpass", "PT Sans
                    Narrow", "Raleway", "Times New Roman".
                lineposition
                    Sets the kind of decoration line(s) with text,
                    such as an "under", "over" or "through" as well
                    as combinations e.g. "under+over", etc.
                shadow
                    Sets the shape and color of the shadow behind
                    text. "auto" places minimal shadow and applies
                    contrast text font color. See
                    https://developer.mozilla.org/en-
                    US/docs/Web/CSS/text-shadow for additional
                    options.
                size

                style
                    Sets whether a font should be styled with a
                    normal or italic face from its family.
                textcase
                    Sets capitalization of text. It can be used to
                    make text appear in all-uppercase or all-
                    lowercase, or with each word capitalized.
                variant
                    Sets the variant of the font.
                weight
                    Sets the weight (or boldness) of the font.

        Returns
        -------
        plotly.graph_objs.layout.legend.Font
        """
        return self["font"]

    @font.setter
    def font(self, val):
        self["font"] = val

    # groupclick
    # ----------
    @property
    def groupclick(self):
        """
        Determines the behavior on legend group item click.
        "toggleitem" toggles the visibility of the individual item
        clicked on the graph. "togglegroup" toggles the visibility of
        all items in the same legendgroup as the item clicked on the
        graph.

        The 'groupclick' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['toggleitem', 'togglegroup']

        Returns
        -------
        Any
        """
        return self["groupclick"]

    @groupclick.setter
    def groupclick(self, val):
        self["groupclick"] = val

    # grouptitlefont
    # --------------
    @property
    def grouptitlefont(self):
        """
        Sets the font for group titles in legend. Defaults to
        `legend.font` with its size increased about 10%.

        The 'grouptitlefont' property is an instance of Grouptitlefont
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.legend.Grouptitlefont`
          - A dict of string/value properties that will be passed
            to the Grouptitlefont constructor

            Supported dict properties:

                color

                family
                    HTML font family - the typeface that will be
                    applied by the web browser. The web browser
                    will only be able to apply a font if it is
                    available on the system which it operates.
                    Provide multiple font families, separated by
                    commas, to indicate the preference in which to
                    apply fonts if they aren't available on the
                    system. The Chart Studio Cloud (at
                    https://chart-studio.plotly.com or on-premise)
                    generates images on a server, where only a
                    select number of fonts are installed and
                    supported. These include "Arial", "Balto",
                    "Courier New", "Droid Sans", "Droid Serif",
                    "Droid Sans Mono", "Gravitas One", "Old
                    Standard TT", "Open Sans", "Overpass", "PT Sans
                    Narrow", "Raleway", "Times New Roman".
                lineposition
                    Sets the kind of decoration line(s) with text,
                    such as an "under", "over" or "through" as well
                    as combinations e.g. "under+over", etc.
                shadow
                    Sets the shape and color of the shadow behind
                    text. "auto" places minimal shadow and applies
                    contrast text font color. See
                    https://developer.mozilla.org/en-
                    US/docs/Web/CSS/text-shadow for additional
                    options.
                size

                style
                    Sets whether a font should be styled with a
                    normal or italic face from its family.
                textcase
                    Sets capitalization of text. It can be used to
                    make text appear in all-uppercase or all-
                    lowercase, or with each word capitalized.
                variant
                    Sets the variant of the font.
                weight
                    Sets the weight (or boldness) of the font.

        Returns
        -------
        plotly.graph_objs.layout.legend.Grouptitlefont
        """
        return self["grouptitlefont"]

    @grouptitlefont.setter
    def grouptitlefont(self, val):
        self["grouptitlefont"] = val

    # indentation
    # -----------
    @property
    def indentation(self):
        """
        Sets the indentation (in px) of the legend entries.

        The 'indentation' property is a number and may be specified as:
          - An int or float in the interval [-15, inf]

        Returns
        -------
        int|float
        """
        return self["indentation"]

    @indentation.setter
    def indentation(self, val):
        self["indentation"] = val

    # itemclick
    # ---------
    @property
    def itemclick(self):
        """
        Determines the behavior on legend item click. "toggle" toggles
        the visibility of the item clicked on the graph. "toggleothers"
        makes the clicked item the sole visible item on the graph.
        False disables legend item click interactions.

        The 'itemclick' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['toggle', 'toggleothers', False]

        Returns
        -------
        Any
        """
        return self["itemclick"]

    @itemclick.setter
    def itemclick(self, val):
        self["itemclick"] = val

    # itemdoubleclick
    # ---------------
    @property
    def itemdoubleclick(self):
        """
        Determines the behavior on legend item double-click. "toggle"
        toggles the visibility of the item clicked on the graph.
        "toggleothers" makes the clicked item the sole visible item on
        the graph. False disables legend item double-click
        interactions.

        The 'itemdoubleclick' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['toggle', 'toggleothers', False]

        Returns
        -------
        Any
        """
        return self["itemdoubleclick"]

    @itemdoubleclick.setter
    def itemdoubleclick(self, val):
        self["itemdoubleclick"] = val

    # itemsizing
    # ----------
    @property
    def itemsizing(self):
        """
        Determines if the legend items symbols scale with their
        corresponding "trace" attributes or remain "constant"
        independent of the symbol size on the graph.

        The 'itemsizing' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['trace', 'constant']

        Returns
        -------
        Any
        """
        return self["itemsizing"]

    @itemsizing.setter
    def itemsizing(self, val):
        self["itemsizing"] = val

    # itemwidth
    # ---------
    @property
    def itemwidth(self):
        """
        Sets the width (in px) of the legend item symbols (the part
        other than the title.text).

        The 'itemwidth' property is a number and may be specified as:
          - An int or float in the interval [30, inf]

        Returns
        -------
        int|float
        """
        return self["itemwidth"]

    @itemwidth.setter
    def itemwidth(self, val):
        self["itemwidth"] = val

    # orientation
    # -----------
    @property
    def orientation(self):
        """
        Sets the orientation of the legend.

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

    # title
    # -----
    @property
    def title(self):
        """
        The 'title' property is an instance of Title
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.legend.Title`
          - A dict of string/value properties that will be passed
            to the Title constructor

            Supported dict properties:

                font
                    Sets this legend's title font. Defaults to
                    `legend.font` with its size increased about
                    20%.
                side
                    Determines the location of legend's title with
                    respect to the legend items. Defaulted to "top"
                    with `orientation` is "h". Defaulted to "left"
                    with `orientation` is "v". The *top left*
                    options could be used to expand top center and
                    top right are for horizontal alignment legend
                    area in both x and y sides.
                text
                    Sets the title of the legend.

        Returns
        -------
        plotly.graph_objs.layout.legend.Title
        """
        return self["title"]

    @title.setter
    def title(self, val):
        self["title"] = val

    # tracegroupgap
    # -------------
    @property
    def tracegroupgap(self):
        """
        Sets the amount of vertical space (in px) between legend
        groups.

        The 'tracegroupgap' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["tracegroupgap"]

    @tracegroupgap.setter
    def tracegroupgap(self, val):
        self["tracegroupgap"] = val

    # traceorder
    # ----------
    @property
    def traceorder(self):
        """
        Determines the order at which the legend items are displayed.
        If "normal", the items are displayed top-to-bottom in the same
        order as the input data. If "reversed", the items are displayed
        in the opposite order as "normal". If "grouped", the items are
        displayed in groups (when a trace `legendgroup` is provided).
        if "grouped+reversed", the items are displayed in the opposite
        order as "grouped".

        The 'traceorder' property is a flaglist and may be specified
        as a string containing:
          - Any combination of ['reversed', 'grouped'] joined with '+' characters
            (e.g. 'reversed+grouped')
            OR exactly one of ['normal'] (e.g. 'normal')

        Returns
        -------
        Any
        """
        return self["traceorder"]

    @traceorder.setter
    def traceorder(self, val):
        self["traceorder"] = val

    # uirevision
    # ----------
    @property
    def uirevision(self):
        """
        Controls persistence of legend-driven changes in trace and pie
        label visibility. Defaults to `layout.uirevision`.

        The 'uirevision' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["uirevision"]

    @uirevision.setter
    def uirevision(self, val):
        self["uirevision"] = val

    # valign
    # ------
    @property
    def valign(self):
        """
        Sets the vertical alignment of the symbols with respect to
        their associated text.

        The 'valign' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['top', 'middle', 'bottom']

        Returns
        -------
        Any
        """
        return self["valign"]

    @valign.setter
    def valign(self, val):
        self["valign"] = val

    # visible
    # -------
    @property
    def visible(self):
        """
        Determines whether or not this legend is visible.

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

    # x
    # -
    @property
    def x(self):
        """
        Sets the x position with respect to `xref` (in normalized
        coordinates) of the legend. When `xref` is "paper", defaults to
        1.02 for vertical legends and defaults to 0 for horizontal
        legends. When `xref` is "container", defaults to 1 for vertical
        legends and defaults to 0 for horizontal legends. Must be
        between 0 and 1 if `xref` is "container". and between "-2" and
        3 if `xref` is "paper".

        The 'x' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["x"]

    @x.setter
    def x(self, val):
        self["x"] = val

    # xanchor
    # -------
    @property
    def xanchor(self):
        """
        Sets the legend's horizontal position anchor. This anchor binds
        the `x` position to the "left", "center" or "right" of the
        legend. Value "auto" anchors legends to the right for `x`
        values greater than or equal to 2/3, anchors legends to the
        left for `x` values less than or equal to 1/3 and anchors
        legends with respect to their center otherwise.

        The 'xanchor' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['auto', 'left', 'center', 'right']

        Returns
        -------
        Any
        """
        return self["xanchor"]

    @xanchor.setter
    def xanchor(self, val):
        self["xanchor"] = val

    # xref
    # ----
    @property
    def xref(self):
        """
        Sets the container `x` refers to. "container" spans the entire
        `width` of the plot. "paper" refers to the width of the
        plotting area only.

        The 'xref' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['container', 'paper']

        Returns
        -------
        Any
        """
        return self["xref"]

    @xref.setter
    def xref(self, val):
        self["xref"] = val

    # y
    # -
    @property
    def y(self):
        """
        Sets the y position with respect to `yref` (in normalized
        coordinates) of the legend. When `yref` is "paper", defaults to
        1 for vertical legends, defaults to "-0.1" for horizontal
        legends on graphs w/o range sliders and defaults to 1.1 for
        horizontal legends on graph with one or multiple range sliders.
        When `yref` is "container", defaults to 1. Must be between 0
        and 1 if `yref` is "container" and between "-2" and 3 if `yref`
        is "paper".

        The 'y' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["y"]

    @y.setter
    def y(self, val):
        self["y"] = val

    # yanchor
    # -------
    @property
    def yanchor(self):
        """
        Sets the legend's vertical position anchor This anchor binds
        the `y` position to the "top", "middle" or "bottom" of the
        legend. Value "auto" anchors legends at their bottom for `y`
        values less than or equal to 1/3, anchors legends to at their
        top for `y` values greater than or equal to 2/3 and anchors
        legends with respect to their middle otherwise.

        The 'yanchor' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['auto', 'top', 'middle', 'bottom']

        Returns
        -------
        Any
        """
        return self["yanchor"]

    @yanchor.setter
    def yanchor(self, val):
        self["yanchor"] = val

    # yref
    # ----
    @property
    def yref(self):
        """
        Sets the container `y` refers to. "container" spans the entire
        `height` of the plot. "paper" refers to the height of the
        plotting area only.

        The 'yref' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['container', 'paper']

        Returns
        -------
        Any
        """
        return self["yref"]

    @yref.setter
    def yref(self, val):
        self["yref"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        bgcolor
            Sets the legend background color. Defaults to
            `layout.paper_bgcolor`.
        bordercolor
            Sets the color of the border enclosing the legend.
        borderwidth
            Sets the width (in px) of the border enclosing the
            legend.
        entrywidth
            Sets the width (in px or fraction) of the legend. Use 0
            to size the entry based on the text width, when
            `entrywidthmode` is set to "pixels".
        entrywidthmode
            Determines what entrywidth means.
        font
            Sets the font used to text the legend items.
        groupclick
            Determines the behavior on legend group item click.
            "toggleitem" toggles the visibility of the individual
            item clicked on the graph. "togglegroup" toggles the
            visibility of all items in the same legendgroup as the
            item clicked on the graph.
        grouptitlefont
            Sets the font for group titles in legend. Defaults to
            `legend.font` with its size increased about 10%.
        indentation
            Sets the indentation (in px) of the legend entries.
        itemclick
            Determines the behavior on legend item click. "toggle"
            toggles the visibility of the item clicked on the
            graph. "toggleothers" makes the clicked item the sole
            visible item on the graph. False disables legend item
            click interactions.
        itemdoubleclick
            Determines the behavior on legend item double-click.
            "toggle" toggles the visibility of the item clicked on
            the graph. "toggleothers" makes the clicked item the
            sole visible item on the graph. False disables legend
            item double-click interactions.
        itemsizing
            Determines if the legend items symbols scale with their
            corresponding "trace" attributes or remain "constant"
            independent of the symbol size on the graph.
        itemwidth
            Sets the width (in px) of the legend item symbols (the
            part other than the title.text).
        orientation
            Sets the orientation of the legend.
        title
            :class:`plotly.graph_objects.layout.legend.Title`
            instance or dict with compatible properties
        tracegroupgap
            Sets the amount of vertical space (in px) between
            legend groups.
        traceorder
            Determines the order at which the legend items are
            displayed. If "normal", the items are displayed top-to-
            bottom in the same order as the input data. If
            "reversed", the items are displayed in the opposite
            order as "normal". If "grouped", the items are
            displayed in groups (when a trace `legendgroup` is
            provided). if "grouped+reversed", the items are
            displayed in the opposite order as "grouped".
        uirevision
            Controls persistence of legend-driven changes in trace
            and pie label visibility. Defaults to
            `layout.uirevision`.
        valign
            Sets the vertical alignment of the symbols with respect
            to their associated text.
        visible
            Determines whether or not this legend is visible.
        x
            Sets the x position with respect to `xref` (in
            normalized coordinates) of the legend. When `xref` is
            "paper", defaults to 1.02 for vertical legends and
            defaults to 0 for horizontal legends. When `xref` is
            "container", defaults to 1 for vertical legends and
            defaults to 0 for horizontal legends. Must be between 0
            and 1 if `xref` is "container". and between "-2" and 3
            if `xref` is "paper".
        xanchor
            Sets the legend's horizontal position anchor. This
            anchor binds the `x` position to the "left", "center"
            or "right" of the legend. Value "auto" anchors legends
            to the right for `x` values greater than or equal to
            2/3, anchors legends to the left for `x` values less
            than or equal to 1/3 and anchors legends with respect
            to their center otherwise.
        xref
            Sets the container `x` refers to. "container" spans the
            entire `width` of the plot. "paper" refers to the width
            of the plotting area only.
        y
            Sets the y position with respect to `yref` (in
            normalized coordinates) of the legend. When `yref` is
            "paper", defaults to 1 for vertical legends, defaults
            to "-0.1" for horizontal legends on graphs w/o range
            sliders and defaults to 1.1 for horizontal legends on
            graph with one or multiple range sliders. When `yref`
            is "container", defaults to 1. Must be between 0 and 1
            if `yref` is "container" and between "-2" and 3 if
            `yref` is "paper".
        yanchor
            Sets the legend's vertical position anchor This anchor
            binds the `y` position to the "top", "middle" or
            "bottom" of the legend. Value "auto" anchors legends at
            their bottom for `y` values less than or equal to 1/3,
            anchors legends to at their top for `y` values greater
            than or equal to 2/3 and anchors legends with respect
            to their middle otherwise.
        yref
            Sets the container `y` refers to. "container" spans the
            entire `height` of the plot. "paper" refers to the
            height of the plotting area only.
        """

    def __init__(
        self,
        arg=None,
        bgcolor=None,
        bordercolor=None,
        borderwidth=None,
        entrywidth=None,
        entrywidthmode=None,
        font=None,
        groupclick=None,
        grouptitlefont=None,
        indentation=None,
        itemclick=None,
        itemdoubleclick=None,
        itemsizing=None,
        itemwidth=None,
        orientation=None,
        title=None,
        tracegroupgap=None,
        traceorder=None,
        uirevision=None,
        valign=None,
        visible=None,
        x=None,
        xanchor=None,
        xref=None,
        y=None,
        yanchor=None,
        yref=None,
        **kwargs,
    ):
        """
        Construct a new Legend object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of :class:`plotly.graph_objs.layout.Legend`
        bgcolor
            Sets the legend background color. Defaults to
            `layout.paper_bgcolor`.
        bordercolor
            Sets the color of the border enclosing the legend.
        borderwidth
            Sets the width (in px) of the border enclosing the
            legend.
        entrywidth
            Sets the width (in px or fraction) of the legend. Use 0
            to size the entry based on the text width, when
            `entrywidthmode` is set to "pixels".
        entrywidthmode
            Determines what entrywidth means.
        font
            Sets the font used to text the legend items.
        groupclick
            Determines the behavior on legend group item click.
            "toggleitem" toggles the visibility of the individual
            item clicked on the graph. "togglegroup" toggles the
            visibility of all items in the same legendgroup as the
            item clicked on the graph.
        grouptitlefont
            Sets the font for group titles in legend. Defaults to
            `legend.font` with its size increased about 10%.
        indentation
            Sets the indentation (in px) of the legend entries.
        itemclick
            Determines the behavior on legend item click. "toggle"
            toggles the visibility of the item clicked on the
            graph. "toggleothers" makes the clicked item the sole
            visible item on the graph. False disables legend item
            click interactions.
        itemdoubleclick
            Determines the behavior on legend item double-click.
            "toggle" toggles the visibility of the item clicked on
            the graph. "toggleothers" makes the clicked item the
            sole visible item on the graph. False disables legend
            item double-click interactions.
        itemsizing
            Determines if the legend items symbols scale with their
            corresponding "trace" attributes or remain "constant"
            independent of the symbol size on the graph.
        itemwidth
            Sets the width (in px) of the legend item symbols (the
            part other than the title.text).
        orientation
            Sets the orientation of the legend.
        title
            :class:`plotly.graph_objects.layout.legend.Title`
            instance or dict with compatible properties
        tracegroupgap
            Sets the amount of vertical space (in px) between
            legend groups.
        traceorder
            Determines the order at which the legend items are
            displayed. If "normal", the items are displayed top-to-
            bottom in the same order as the input data. If
            "reversed", the items are displayed in the opposite
            order as "normal". If "grouped", the items are
            displayed in groups (when a trace `legendgroup` is
            provided). if "grouped+reversed", the items are
            displayed in the opposite order as "grouped".
        uirevision
            Controls persistence of legend-driven changes in trace
            and pie label visibility. Defaults to
            `layout.uirevision`.
        valign
            Sets the vertical alignment of the symbols with respect
            to their associated text.
        visible
            Determines whether or not this legend is visible.
        x
            Sets the x position with respect to `xref` (in
            normalized coordinates) of the legend. When `xref` is
            "paper", defaults to 1.02 for vertical legends and
            defaults to 0 for horizontal legends. When `xref` is
            "container", defaults to 1 for vertical legends and
            defaults to 0 for horizontal legends. Must be between 0
            and 1 if `xref` is "container". and between "-2" and 3
            if `xref` is "paper".
        xanchor
            Sets the legend's horizontal position anchor. This
            anchor binds the `x` position to the "left", "center"
            or "right" of the legend. Value "auto" anchors legends
            to the right for `x` values greater than or equal to
            2/3, anchors legends to the left for `x` values less
            than or equal to 1/3 and anchors legends with respect
            to their center otherwise.
        xref
            Sets the container `x` refers to. "container" spans the
            entire `width` of the plot. "paper" refers to the width
            of the plotting area only.
        y
            Sets the y position with respect to `yref` (in
            normalized coordinates) of the legend. When `yref` is
            "paper", defaults to 1 for vertical legends, defaults
            to "-0.1" for horizontal legends on graphs w/o range
            sliders and defaults to 1.1 for horizontal legends on
            graph with one or multiple range sliders. When `yref`
            is "container", defaults to 1. Must be between 0 and 1
            if `yref` is "container" and between "-2" and 3 if
            `yref` is "paper".
        yanchor
            Sets the legend's vertical position anchor This anchor
            binds the `y` position to the "top", "middle" or
            "bottom" of the legend. Value "auto" anchors legends at
            their bottom for `y` values less than or equal to 1/3,
            anchors legends to at their top for `y` values greater
            than or equal to 2/3 and anchors legends with respect
            to their middle otherwise.
        yref
            Sets the container `y` refers to. "container" spans the
            entire `height` of the plot. "paper" refers to the
            height of the plotting area only.

        Returns
        -------
        Legend
        """
        super(Legend, self).__init__("legend")

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
The first argument to the plotly.graph_objs.layout.Legend
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.Legend`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("bgcolor", None)
        _v = bgcolor if bgcolor is not None else _v
        if _v is not None:
            self["bgcolor"] = _v
        _v = arg.pop("bordercolor", None)
        _v = bordercolor if bordercolor is not None else _v
        if _v is not None:
            self["bordercolor"] = _v
        _v = arg.pop("borderwidth", None)
        _v = borderwidth if borderwidth is not None else _v
        if _v is not None:
            self["borderwidth"] = _v
        _v = arg.pop("entrywidth", None)
        _v = entrywidth if entrywidth is not None else _v
        if _v is not None:
            self["entrywidth"] = _v
        _v = arg.pop("entrywidthmode", None)
        _v = entrywidthmode if entrywidthmode is not None else _v
        if _v is not None:
            self["entrywidthmode"] = _v
        _v = arg.pop("font", None)
        _v = font if font is not None else _v
        if _v is not None:
            self["font"] = _v
        _v = arg.pop("groupclick", None)
        _v = groupclick if groupclick is not None else _v
        if _v is not None:
            self["groupclick"] = _v
        _v = arg.pop("grouptitlefont", None)
        _v = grouptitlefont if grouptitlefont is not None else _v
        if _v is not None:
            self["grouptitlefont"] = _v
        _v = arg.pop("indentation", None)
        _v = indentation if indentation is not None else _v
        if _v is not None:
            self["indentation"] = _v
        _v = arg.pop("itemclick", None)
        _v = itemclick if itemclick is not None else _v
        if _v is not None:
            self["itemclick"] = _v
        _v = arg.pop("itemdoubleclick", None)
        _v = itemdoubleclick if itemdoubleclick is not None else _v
        if _v is not None:
            self["itemdoubleclick"] = _v
        _v = arg.pop("itemsizing", None)
        _v = itemsizing if itemsizing is not None else _v
        if _v is not None:
            self["itemsizing"] = _v
        _v = arg.pop("itemwidth", None)
        _v = itemwidth if itemwidth is not None else _v
        if _v is not None:
            self["itemwidth"] = _v
        _v = arg.pop("orientation", None)
        _v = orientation if orientation is not None else _v
        if _v is not None:
            self["orientation"] = _v
        _v = arg.pop("title", None)
        _v = title if title is not None else _v
        if _v is not None:
            self["title"] = _v
        _v = arg.pop("tracegroupgap", None)
        _v = tracegroupgap if tracegroupgap is not None else _v
        if _v is not None:
            self["tracegroupgap"] = _v
        _v = arg.pop("traceorder", None)
        _v = traceorder if traceorder is not None else _v
        if _v is not None:
            self["traceorder"] = _v
        _v = arg.pop("uirevision", None)
        _v = uirevision if uirevision is not None else _v
        if _v is not None:
            self["uirevision"] = _v
        _v = arg.pop("valign", None)
        _v = valign if valign is not None else _v
        if _v is not None:
            self["valign"] = _v
        _v = arg.pop("visible", None)
        _v = visible if visible is not None else _v
        if _v is not None:
            self["visible"] = _v
        _v = arg.pop("x", None)
        _v = x if x is not None else _v
        if _v is not None:
            self["x"] = _v
        _v = arg.pop("xanchor", None)
        _v = xanchor if xanchor is not None else _v
        if _v is not None:
            self["xanchor"] = _v
        _v = arg.pop("xref", None)
        _v = xref if xref is not None else _v
        if _v is not None:
            self["xref"] = _v
        _v = arg.pop("y", None)
        _v = y if y is not None else _v
        if _v is not None:
            self["y"] = _v
        _v = arg.pop("yanchor", None)
        _v = yanchor if yanchor is not None else _v
        if _v is not None:
            self["yanchor"] = _v
        _v = arg.pop("yref", None)
        _v = yref if yref is not None else _v
        if _v is not None:
            self["yref"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
