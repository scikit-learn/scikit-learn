#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Legend(_BaseLayoutHierarchyType):
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
        "maxheight",
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
    def bordercolor(self):
        """
        Sets the color of the border enclosing the legend.

        The 'bordercolor' property is a color and may be specified as:
          - A hex string (e.g. '#ff0000')
          - An rgb/rgba string (e.g. 'rgb(255,0,0)')
          - An hsl/hsla string (e.g. 'hsl(0,100%,50%)')
          - An hsv/hsva string (e.g. 'hsv(0,100%,100%)')
          - A named CSS color: see https://plotly.com/python/css-colors/ for a list

        Returns
        -------
        str
        """
        return self["bordercolor"]

    @bordercolor.setter
    def bordercolor(self, val):
        self["bordercolor"] = val

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

    @property
    def font(self):
        """
        Sets the font used to text the legend items.

        The 'font' property is an instance of Font
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.legend.Font`
          - A dict of string/value properties that will be passed
            to the Font constructor

        Returns
        -------
        plotly.graph_objs.layout.legend.Font
        """
        return self["font"]

    @font.setter
    def font(self, val):
        self["font"] = val

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

        Returns
        -------
        plotly.graph_objs.layout.legend.Grouptitlefont
        """
        return self["grouptitlefont"]

    @grouptitlefont.setter
    def grouptitlefont(self, val):
        self["grouptitlefont"] = val

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

    @property
    def maxheight(self):
        """
        Sets the max height (in px) of the legend, or max height ratio
        (reference height * ratio) if less than or equal to 1. Default
        value is: 0.5 for horizontal legends; 1 for vertical legends.
        The minimum allowed height is 30px. For a ratio of 0.5, the
        legend will take up to 50% of the reference height before
        displaying a scrollbar. The reference height is the full layout
        height with the following exception: vertically oriented
        legends with a `yref` of `"paper", located to the side of the
        plot. In this case, the reference height is the plot height.

        The 'maxheight' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["maxheight"]

    @maxheight.setter
    def maxheight(self, val):
        self["maxheight"] = val

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

    @property
    def title(self):
        """
        The 'title' property is an instance of Title
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.legend.Title`
          - A dict of string/value properties that will be passed
            to the Title constructor

        Returns
        -------
        plotly.graph_objs.layout.legend.Title
        """
        return self["title"]

    @title.setter
    def title(self, val):
        self["title"] = val

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

    @property
    def yanchor(self):
        """
        Sets the legend's vertical position anchor. This anchor binds
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
        maxheight
            Sets the max height (in px) of the legend, or max
            height ratio (reference height * ratio) if less than or
            equal to 1. Default value is: 0.5 for horizontal
            legends; 1 for vertical legends. The minimum allowed
            height is 30px. For a ratio of 0.5, the legend will
            take up to 50% of the reference height before
            displaying a scrollbar. The reference height is the
            full layout height with the following exception:
            vertically oriented legends with a `yref` of `"paper",
            located to the side of the plot. In this case, the
            reference height is the plot height.
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
            Sets the legend's vertical position anchor. This anchor
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
        maxheight=None,
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
        maxheight
            Sets the max height (in px) of the legend, or max
            height ratio (reference height * ratio) if less than or
            equal to 1. Default value is: 0.5 for horizontal
            legends; 1 for vertical legends. The minimum allowed
            height is 30px. For a ratio of 0.5, the legend will
            take up to 50% of the reference height before
            displaying a scrollbar. The reference height is the
            full layout height with the following exception:
            vertically oriented legends with a `yref` of `"paper",
            located to the side of the plot. In this case, the
            reference height is the plot height.
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
            Sets the legend's vertical position anchor. This anchor
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
        super().__init__("legend")
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
The first argument to the plotly.graph_objs.layout.Legend
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.Legend`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("bgcolor", arg, bgcolor)
        self._set_property("bordercolor", arg, bordercolor)
        self._set_property("borderwidth", arg, borderwidth)
        self._set_property("entrywidth", arg, entrywidth)
        self._set_property("entrywidthmode", arg, entrywidthmode)
        self._set_property("font", arg, font)
        self._set_property("groupclick", arg, groupclick)
        self._set_property("grouptitlefont", arg, grouptitlefont)
        self._set_property("indentation", arg, indentation)
        self._set_property("itemclick", arg, itemclick)
        self._set_property("itemdoubleclick", arg, itemdoubleclick)
        self._set_property("itemsizing", arg, itemsizing)
        self._set_property("itemwidth", arg, itemwidth)
        self._set_property("maxheight", arg, maxheight)
        self._set_property("orientation", arg, orientation)
        self._set_property("title", arg, title)
        self._set_property("tracegroupgap", arg, tracegroupgap)
        self._set_property("traceorder", arg, traceorder)
        self._set_property("uirevision", arg, uirevision)
        self._set_property("valign", arg, valign)
        self._set_property("visible", arg, visible)
        self._set_property("x", arg, x)
        self._set_property("xanchor", arg, xanchor)
        self._set_property("xref", arg, xref)
        self._set_property("y", arg, y)
        self._set_property("yanchor", arg, yanchor)
        self._set_property("yref", arg, yref)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
