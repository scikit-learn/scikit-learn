from plotly.basedatatypes import BaseLayoutType as _BaseLayoutType
import copy as _copy


class Layout(_BaseLayoutType):

    _subplotid_prop_names = [
        "coloraxis",
        "geo",
        "legend",
        "map",
        "mapbox",
        "polar",
        "scene",
        "smith",
        "ternary",
        "xaxis",
        "yaxis",
    ]

    import re

    _subplotid_prop_re = re.compile("^(" + "|".join(_subplotid_prop_names) + r")(\d+)$")

    @property
    def _subplotid_validators(self):
        """
        dict of validator classes for each subplot type

        Returns
        -------
        dict
        """
        from plotly.validators.layout import (
            ColoraxisValidator,
            GeoValidator,
            LegendValidator,
            MapValidator,
            MapboxValidator,
            PolarValidator,
            SceneValidator,
            SmithValidator,
            TernaryValidator,
            XaxisValidator,
            YaxisValidator,
        )

        return {
            "coloraxis": ColoraxisValidator,
            "geo": GeoValidator,
            "legend": LegendValidator,
            "map": MapValidator,
            "mapbox": MapboxValidator,
            "polar": PolarValidator,
            "scene": SceneValidator,
            "smith": SmithValidator,
            "ternary": TernaryValidator,
            "xaxis": XaxisValidator,
            "yaxis": YaxisValidator,
        }

    def _subplot_re_match(self, prop):
        return self._subplotid_prop_re.match(prop)

    # class properties
    # --------------------
    _parent_path_str = ""
    _path_str = "layout"
    _valid_props = {
        "activeselection",
        "activeshape",
        "annotationdefaults",
        "annotations",
        "autosize",
        "autotypenumbers",
        "barcornerradius",
        "bargap",
        "bargroupgap",
        "barmode",
        "barnorm",
        "boxgap",
        "boxgroupgap",
        "boxmode",
        "calendar",
        "clickmode",
        "coloraxis",
        "colorscale",
        "colorway",
        "computed",
        "datarevision",
        "dragmode",
        "editrevision",
        "extendfunnelareacolors",
        "extendiciclecolors",
        "extendpiecolors",
        "extendsunburstcolors",
        "extendtreemapcolors",
        "font",
        "funnelareacolorway",
        "funnelgap",
        "funnelgroupgap",
        "funnelmode",
        "geo",
        "grid",
        "height",
        "hiddenlabels",
        "hiddenlabelssrc",
        "hidesources",
        "hoverdistance",
        "hoverlabel",
        "hovermode",
        "hoversubplots",
        "iciclecolorway",
        "imagedefaults",
        "images",
        "legend",
        "map",
        "mapbox",
        "margin",
        "meta",
        "metasrc",
        "minreducedheight",
        "minreducedwidth",
        "modebar",
        "newselection",
        "newshape",
        "paper_bgcolor",
        "piecolorway",
        "plot_bgcolor",
        "polar",
        "scattergap",
        "scattermode",
        "scene",
        "selectdirection",
        "selectiondefaults",
        "selectionrevision",
        "selections",
        "separators",
        "shapedefaults",
        "shapes",
        "showlegend",
        "sliderdefaults",
        "sliders",
        "smith",
        "spikedistance",
        "sunburstcolorway",
        "template",
        "ternary",
        "title",
        "titlefont",
        "transition",
        "treemapcolorway",
        "uirevision",
        "uniformtext",
        "updatemenudefaults",
        "updatemenus",
        "violingap",
        "violingroupgap",
        "violinmode",
        "waterfallgap",
        "waterfallgroupgap",
        "waterfallmode",
        "width",
        "xaxis",
        "yaxis",
    }

    # activeselection
    # ---------------
    @property
    def activeselection(self):
        """
        The 'activeselection' property is an instance of Activeselection
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Activeselection`
          - A dict of string/value properties that will be passed
            to the Activeselection constructor

            Supported dict properties:

                fillcolor
                    Sets the color filling the active selection'
                    interior.
                opacity
                    Sets the opacity of the active selection.

        Returns
        -------
        plotly.graph_objs.layout.Activeselection
        """
        return self["activeselection"]

    @activeselection.setter
    def activeselection(self, val):
        self["activeselection"] = val

    # activeshape
    # -----------
    @property
    def activeshape(self):
        """
        The 'activeshape' property is an instance of Activeshape
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Activeshape`
          - A dict of string/value properties that will be passed
            to the Activeshape constructor

            Supported dict properties:

                fillcolor
                    Sets the color filling the active shape'
                    interior.
                opacity
                    Sets the opacity of the active shape.

        Returns
        -------
        plotly.graph_objs.layout.Activeshape
        """
        return self["activeshape"]

    @activeshape.setter
    def activeshape(self, val):
        self["activeshape"] = val

    # annotations
    # -----------
    @property
    def annotations(self):
        """
        The 'annotations' property is a tuple of instances of
        Annotation that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.Annotation
          - A list or tuple of dicts of string/value properties that
            will be passed to the Annotation constructor

            Supported dict properties:

                align
                    Sets the horizontal alignment of the `text`
                    within the box. Has an effect only if `text`
                    spans two or more lines (i.e. `text` contains
                    one or more <br> HTML tags) or if an explicit
                    width is set to override the text width.
                arrowcolor
                    Sets the color of the annotation arrow.
                arrowhead
                    Sets the end annotation arrow head style.
                arrowside
                    Sets the annotation arrow head position.
                arrowsize
                    Sets the size of the end annotation arrow head,
                    relative to `arrowwidth`. A value of 1
                    (default) gives a head about 3x as wide as the
                    line.
                arrowwidth
                    Sets the width (in px) of annotation arrow
                    line.
                ax
                    Sets the x component of the arrow tail about
                    the arrow head. If `axref` is `pixel`, a
                    positive (negative) component corresponds to an
                    arrow pointing from right to left (left to
                    right). If `axref` is not `pixel` and is
                    exactly the same as `xref`, this is an absolute
                    value on that axis, like `x`, specified in the
                    same coordinates as `xref`.
                axref
                    Indicates in what coordinates the tail of the
                    annotation (ax,ay) is specified. If set to a x
                    axis id (e.g. "x" or "x2"), the `x` position
                    refers to a x coordinate. If set to "paper",
                    the `x` position refers to the distance from
                    the left of the plotting area in normalized
                    coordinates where 0 (1) corresponds to the left
                    (right). If set to a x axis ID followed by
                    "domain" (separated by a space), the position
                    behaves like for "paper", but refers to the
                    distance in fractions of the domain length from
                    the left of the domain of that axis: e.g., *x2
                    domain* refers to the domain of the second x
                    axis and a x position of 0.5 refers to the
                    point between the left and the right of the
                    domain of the second x axis. In order for
                    absolute positioning of the arrow to work,
                    "axref" must be exactly the same as "xref",
                    otherwise "axref" will revert to "pixel"
                    (explained next). For relative positioning,
                    "axref" can be set to "pixel", in which case
                    the "ax" value is specified in pixels relative
                    to "x". Absolute positioning is useful for
                    trendline annotations which should continue to
                    indicate the correct trend when zoomed.
                    Relative positioning is useful for specifying
                    the text offset for an annotated point.
                ay
                    Sets the y component of the arrow tail about
                    the arrow head. If `ayref` is `pixel`, a
                    positive (negative) component corresponds to an
                    arrow pointing from bottom to top (top to
                    bottom). If `ayref` is not `pixel` and is
                    exactly the same as `yref`, this is an absolute
                    value on that axis, like `y`, specified in the
                    same coordinates as `yref`.
                ayref
                    Indicates in what coordinates the tail of the
                    annotation (ax,ay) is specified. If set to a y
                    axis id (e.g. "y" or "y2"), the `y` position
                    refers to a y coordinate. If set to "paper",
                    the `y` position refers to the distance from
                    the bottom of the plotting area in normalized
                    coordinates where 0 (1) corresponds to the
                    bottom (top). If set to a y axis ID followed by
                    "domain" (separated by a space), the position
                    behaves like for "paper", but refers to the
                    distance in fractions of the domain length from
                    the bottom of the domain of that axis: e.g.,
                    *y2 domain* refers to the domain of the second
                    y  axis and a y position of 0.5 refers to the
                    point between the bottom and the top of the
                    domain of the second y axis. In order for
                    absolute positioning of the arrow to work,
                    "ayref" must be exactly the same as "yref",
                    otherwise "ayref" will revert to "pixel"
                    (explained next). For relative positioning,
                    "ayref" can be set to "pixel", in which case
                    the "ay" value is specified in pixels relative
                    to "y". Absolute positioning is useful for
                    trendline annotations which should continue to
                    indicate the correct trend when zoomed.
                    Relative positioning is useful for specifying
                    the text offset for an annotated point.
                bgcolor
                    Sets the background color of the annotation.
                bordercolor
                    Sets the color of the border enclosing the
                    annotation `text`.
                borderpad
                    Sets the padding (in px) between the `text` and
                    the enclosing border.
                borderwidth
                    Sets the width (in px) of the border enclosing
                    the annotation `text`.
                captureevents
                    Determines whether the annotation text box
                    captures mouse move and click events, or allows
                    those events to pass through to data points in
                    the plot that may be behind the annotation. By
                    default `captureevents` is False unless
                    `hovertext` is provided. If you use the event
                    `plotly_clickannotation` without `hovertext`
                    you must explicitly enable `captureevents`.
                clicktoshow
                    Makes this annotation respond to clicks on the
                    plot. If you click a data point that exactly
                    matches the `x` and `y` values of this
                    annotation, and it is hidden (visible: false),
                    it will appear. In "onoff" mode, you must click
                    the same point again to make it disappear, so
                    if you click multiple points, you can show
                    multiple annotations. In "onout" mode, a click
                    anywhere else in the plot (on another data
                    point or not) will hide this annotation. If you
                    need to show/hide this annotation in response
                    to different `x` or `y` values, you can set
                    `xclick` and/or `yclick`. This is useful for
                    example to label the side of a bar. To label
                    markers though, `standoff` is preferred over
                    `xclick` and `yclick`.
                font
                    Sets the annotation text font.
                height
                    Sets an explicit height for the text box. null
                    (default) lets the text set the box height.
                    Taller text will be clipped.
                hoverlabel
                    :class:`plotly.graph_objects.layout.annotation.
                    Hoverlabel` instance or dict with compatible
                    properties
                hovertext
                    Sets text to appear when hovering over this
                    annotation. If omitted or blank, no hover label
                    will appear.
                name
                    When used in a template, named items are
                    created in the output figure in addition to any
                    items the figure already has in this array. You
                    can modify these items in the output figure by
                    making your own item with `templateitemname`
                    matching this `name` alongside your
                    modifications (including `visible: false` or
                    `enabled: false` to hide it). Has no effect
                    outside of a template.
                opacity
                    Sets the opacity of the annotation (text +
                    arrow).
                showarrow
                    Determines whether or not the annotation is
                    drawn with an arrow. If True, `text` is placed
                    near the arrow's tail. If False, `text` lines
                    up with the `x` and `y` provided.
                standoff
                    Sets a distance, in pixels, to move the end
                    arrowhead away from the position it is pointing
                    at, for example to point at the edge of a
                    marker independent of zoom. Note that this
                    shortens the arrow from the `ax` / `ay` vector,
                    in contrast to `xshift` / `yshift` which moves
                    everything by this amount.
                startarrowhead
                    Sets the start annotation arrow head style.
                startarrowsize
                    Sets the size of the start annotation arrow
                    head, relative to `arrowwidth`. A value of 1
                    (default) gives a head about 3x as wide as the
                    line.
                startstandoff
                    Sets a distance, in pixels, to move the start
                    arrowhead away from the position it is pointing
                    at, for example to point at the edge of a
                    marker independent of zoom. Note that this
                    shortens the arrow from the `ax` / `ay` vector,
                    in contrast to `xshift` / `yshift` which moves
                    everything by this amount.
                templateitemname
                    Used to refer to a named item in this array in
                    the template. Named items from the template
                    will be created even without a matching item in
                    the input figure, but you can modify one by
                    making an item with `templateitemname` matching
                    its `name`, alongside your modifications
                    (including `visible: false` or `enabled: false`
                    to hide it). If there is no template or no
                    matching item, this item will be hidden unless
                    you explicitly show it with `visible: true`.
                text
                    Sets the text associated with this annotation.
                    Plotly uses a subset of HTML tags to do things
                    like newline (<br>), bold (<b></b>), italics
                    (<i></i>), hyperlinks (<a href='...'></a>).
                    Tags <em>, <sup>, <sub>, <s>, <u> <span> are
                    also supported.
                textangle
                    Sets the angle at which the `text` is drawn
                    with respect to the horizontal.
                valign
                    Sets the vertical alignment of the `text`
                    within the box. Has an effect only if an
                    explicit height is set to override the text
                    height.
                visible
                    Determines whether or not this annotation is
                    visible.
                width
                    Sets an explicit width for the text box. null
                    (default) lets the text set the box width.
                    Wider text will be clipped. There is no
                    automatic wrapping; use <br> to start a new
                    line.
                x
                    Sets the annotation's x position. If the axis
                    `type` is "log", then you must take the log of
                    your desired range. If the axis `type` is
                    "date", it should be date strings, like date
                    data, though Date objects and unix milliseconds
                    will be accepted and converted to strings. If
                    the axis `type` is "category", it should be
                    numbers, using the scale where each category is
                    assigned a serial number from zero in the order
                    it appears.
                xanchor
                    Sets the text box's horizontal position anchor
                    This anchor binds the `x` position to the
                    "left", "center" or "right" of the annotation.
                    For example, if `x` is set to 1, `xref` to
                    "paper" and `xanchor` to "right" then the
                    right-most portion of the annotation lines up
                    with the right-most edge of the plotting area.
                    If "auto", the anchor is equivalent to "center"
                    for data-referenced annotations or if there is
                    an arrow, whereas for paper-referenced with no
                    arrow, the anchor picked corresponds to the
                    closest side.
                xclick
                    Toggle this annotation when clicking a data
                    point whose `x` value is `xclick` rather than
                    the annotation's `x` value.
                xref
                    Sets the annotation's x coordinate axis. If set
                    to a x axis id (e.g. "x" or "x2"), the `x`
                    position refers to a x coordinate. If set to
                    "paper", the `x` position refers to the
                    distance from the left of the plotting area in
                    normalized coordinates where 0 (1) corresponds
                    to the left (right). If set to a x axis ID
                    followed by "domain" (separated by a space),
                    the position behaves like for "paper", but
                    refers to the distance in fractions of the
                    domain length from the left of the domain of
                    that axis: e.g., *x2 domain* refers to the
                    domain of the second x  axis and a x position
                    of 0.5 refers to the point between the left and
                    the right of the domain of the second x axis.
                xshift
                    Shifts the position of the whole annotation and
                    arrow to the right (positive) or left
                    (negative) by this many pixels.
                y
                    Sets the annotation's y position. If the axis
                    `type` is "log", then you must take the log of
                    your desired range. If the axis `type` is
                    "date", it should be date strings, like date
                    data, though Date objects and unix milliseconds
                    will be accepted and converted to strings. If
                    the axis `type` is "category", it should be
                    numbers, using the scale where each category is
                    assigned a serial number from zero in the order
                    it appears.
                yanchor
                    Sets the text box's vertical position anchor
                    This anchor binds the `y` position to the
                    "top", "middle" or "bottom" of the annotation.
                    For example, if `y` is set to 1, `yref` to
                    "paper" and `yanchor` to "top" then the top-
                    most portion of the annotation lines up with
                    the top-most edge of the plotting area. If
                    "auto", the anchor is equivalent to "middle"
                    for data-referenced annotations or if there is
                    an arrow, whereas for paper-referenced with no
                    arrow, the anchor picked corresponds to the
                    closest side.
                yclick
                    Toggle this annotation when clicking a data
                    point whose `y` value is `yclick` rather than
                    the annotation's `y` value.
                yref
                    Sets the annotation's y coordinate axis. If set
                    to a y axis id (e.g. "y" or "y2"), the `y`
                    position refers to a y coordinate. If set to
                    "paper", the `y` position refers to the
                    distance from the bottom of the plotting area
                    in normalized coordinates where 0 (1)
                    corresponds to the bottom (top). If set to a y
                    axis ID followed by "domain" (separated by a
                    space), the position behaves like for "paper",
                    but refers to the distance in fractions of the
                    domain length from the bottom of the domain of
                    that axis: e.g., *y2 domain* refers to the
                    domain of the second y  axis and a y position
                    of 0.5 refers to the point between the bottom
                    and the top of the domain of the second y axis.
                yshift
                    Shifts the position of the whole annotation and
                    arrow up (positive) or down (negative) by this
                    many pixels.

        Returns
        -------
        tuple[plotly.graph_objs.layout.Annotation]
        """
        return self["annotations"]

    @annotations.setter
    def annotations(self, val):
        self["annotations"] = val

    # annotationdefaults
    # ------------------
    @property
    def annotationdefaults(self):
        """
        When used in a template (as
        layout.template.layout.annotationdefaults), sets the default
        property values to use for elements of layout.annotations

        The 'annotationdefaults' property is an instance of Annotation
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Annotation`
          - A dict of string/value properties that will be passed
            to the Annotation constructor

            Supported dict properties:

        Returns
        -------
        plotly.graph_objs.layout.Annotation
        """
        return self["annotationdefaults"]

    @annotationdefaults.setter
    def annotationdefaults(self, val):
        self["annotationdefaults"] = val

    # autosize
    # --------
    @property
    def autosize(self):
        """
        Determines whether or not a layout width or height that has
        been left undefined by the user is initialized on each
        relayout. Note that, regardless of this attribute, an undefined
        layout width or height is always initialized on the first call
        to plot.

        The 'autosize' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["autosize"]

    @autosize.setter
    def autosize(self, val):
        self["autosize"] = val

    # autotypenumbers
    # ---------------
    @property
    def autotypenumbers(self):
        """
        Using "strict" a numeric string in trace data is not converted
        to a number. Using *convert types* a numeric string in trace
        data may be treated as a number during automatic axis `type`
        detection. This is the default value; however it could be
        overridden for individual axes.

        The 'autotypenumbers' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['convert types', 'strict']

        Returns
        -------
        Any
        """
        return self["autotypenumbers"]

    @autotypenumbers.setter
    def autotypenumbers(self, val):
        self["autotypenumbers"] = val

    # barcornerradius
    # ---------------
    @property
    def barcornerradius(self):
        """
        Sets the rounding of bar corners. May be an integer number of
        pixels, or a percentage of bar width (as a string ending in %).

        The 'barcornerradius' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["barcornerradius"]

    @barcornerradius.setter
    def barcornerradius(self, val):
        self["barcornerradius"] = val

    # bargap
    # ------
    @property
    def bargap(self):
        """
        Sets the gap (in plot fraction) between bars of adjacent
        location coordinates.

        The 'bargap' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["bargap"]

    @bargap.setter
    def bargap(self, val):
        self["bargap"] = val

    # bargroupgap
    # -----------
    @property
    def bargroupgap(self):
        """
        Sets the gap (in plot fraction) between bars of the same
        location coordinate.

        The 'bargroupgap' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["bargroupgap"]

    @bargroupgap.setter
    def bargroupgap(self, val):
        self["bargroupgap"] = val

    # barmode
    # -------
    @property
    def barmode(self):
        """
        Determines how bars at the same location coordinate are
        displayed on the graph. With "stack", the bars are stacked on
        top of one another With "relative", the bars are stacked on top
        of one another, with negative values below the axis, positive
        values above With "group", the bars are plotted next to one
        another centered around the shared location. With "overlay",
        the bars are plotted over one another, you might need to reduce
        "opacity" to see multiple bars.

        The 'barmode' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['stack', 'group', 'overlay', 'relative']

        Returns
        -------
        Any
        """
        return self["barmode"]

    @barmode.setter
    def barmode(self, val):
        self["barmode"] = val

    # barnorm
    # -------
    @property
    def barnorm(self):
        """
        Sets the normalization for bar traces on the graph. With
        "fraction", the value of each bar is divided by the sum of all
        values at that location coordinate. "percent" is the same but
        multiplied by 100 to show percentages.

        The 'barnorm' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['', 'fraction', 'percent']

        Returns
        -------
        Any
        """
        return self["barnorm"]

    @barnorm.setter
    def barnorm(self, val):
        self["barnorm"] = val

    # boxgap
    # ------
    @property
    def boxgap(self):
        """
        Sets the gap (in plot fraction) between boxes of adjacent
        location coordinates. Has no effect on traces that have "width"
        set.

        The 'boxgap' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["boxgap"]

    @boxgap.setter
    def boxgap(self, val):
        self["boxgap"] = val

    # boxgroupgap
    # -----------
    @property
    def boxgroupgap(self):
        """
        Sets the gap (in plot fraction) between boxes of the same
        location coordinate. Has no effect on traces that have "width"
        set.

        The 'boxgroupgap' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["boxgroupgap"]

    @boxgroupgap.setter
    def boxgroupgap(self, val):
        self["boxgroupgap"] = val

    # boxmode
    # -------
    @property
    def boxmode(self):
        """
        Determines how boxes at the same location coordinate are
        displayed on the graph. If "group", the boxes are plotted next
        to one another centered around the shared location. If
        "overlay", the boxes are plotted over one another, you might
        need to set "opacity" to see them multiple boxes. Has no effect
        on traces that have "width" set.

        The 'boxmode' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['group', 'overlay']

        Returns
        -------
        Any
        """
        return self["boxmode"]

    @boxmode.setter
    def boxmode(self, val):
        self["boxmode"] = val

    # calendar
    # --------
    @property
    def calendar(self):
        """
        Sets the default calendar system to use for interpreting and
        displaying dates throughout the plot.

        The 'calendar' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['chinese', 'coptic', 'discworld', 'ethiopian',
                'gregorian', 'hebrew', 'islamic', 'jalali', 'julian',
                'mayan', 'nanakshahi', 'nepali', 'persian', 'taiwan',
                'thai', 'ummalqura']

        Returns
        -------
        Any
        """
        return self["calendar"]

    @calendar.setter
    def calendar(self, val):
        self["calendar"] = val

    # clickmode
    # ---------
    @property
    def clickmode(self):
        """
        Determines the mode of single click interactions. "event" is
        the default value and emits the `plotly_click` event. In
        addition this mode emits the `plotly_selected` event in drag
        modes "lasso" and "select", but with no event data attached
        (kept for compatibility reasons). The "select" flag enables
        selecting single data points via click. This mode also supports
        persistent selections, meaning that pressing Shift while
        clicking, adds to / subtracts from an existing selection.
        "select" with `hovermode`: "x" can be confusing, consider
        explicitly setting `hovermode`: "closest" when using this
        feature. Selection events are sent accordingly as long as
        "event" flag is set as well. When the "event" flag is missing,
        `plotly_click` and `plotly_selected` events are not fired.

        The 'clickmode' property is a flaglist and may be specified
        as a string containing:
          - Any combination of ['event', 'select'] joined with '+' characters
            (e.g. 'event+select')
            OR exactly one of ['none'] (e.g. 'none')

        Returns
        -------
        Any
        """
        return self["clickmode"]

    @clickmode.setter
    def clickmode(self, val):
        self["clickmode"] = val

    # coloraxis
    # ---------
    @property
    def coloraxis(self):
        """
        The 'coloraxis' property is an instance of Coloraxis
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Coloraxis`
          - A dict of string/value properties that will be passed
            to the Coloraxis constructor

            Supported dict properties:

                autocolorscale
                    Determines whether the colorscale is a default
                    palette (`autocolorscale: true`) or the palette
                    determined by `colorscale`. In case
                    `colorscale` is unspecified or `autocolorscale`
                    is true, the default palette will be chosen
                    according to whether numbers in the `color`
                    array are all positive, all negative or mixed.
                cauto
                    Determines whether or not the color domain is
                    computed with respect to the input data (here
                    corresponding trace color array(s)) or the
                    bounds set in `cmin` and `cmax` Defaults to
                    `false` when `cmin` and `cmax` are set by the
                    user.
                cmax
                    Sets the upper bound of the color domain. Value
                    should have the same units as corresponding
                    trace color array(s) and if set, `cmin` must be
                    set as well.
                cmid
                    Sets the mid-point of the color domain by
                    scaling `cmin` and/or `cmax` to be equidistant
                    to this point. Value should have the same units
                    as corresponding trace color array(s). Has no
                    effect when `cauto` is `false`.
                cmin
                    Sets the lower bound of the color domain. Value
                    should have the same units as corresponding
                    trace color array(s) and if set, `cmax` must be
                    set as well.
                colorbar
                    :class:`plotly.graph_objects.layout.coloraxis.C
                    olorBar` instance or dict with compatible
                    properties
                colorscale
                    Sets the colorscale. The colorscale must be an
                    array containing arrays mapping a normalized
                    value to an rgb, rgba, hex, hsl, hsv, or named
                    color string. At minimum, a mapping for the
                    lowest (0) and highest (1) values are required.
                    For example, `[[0, 'rgb(0,0,255)'], [1,
                    'rgb(255,0,0)']]`. To control the bounds of the
                    colorscale in color space, use `cmin` and
                    `cmax`. Alternatively, `colorscale` may be a
                    palette name string of the following list: Blac
                    kbody,Bluered,Blues,Cividis,Earth,Electric,Gree
                    ns,Greys,Hot,Jet,Picnic,Portland,Rainbow,RdBu,R
                    eds,Viridis,YlGnBu,YlOrRd.
                reversescale
                    Reverses the color mapping if true. If true,
                    `cmin` will correspond to the last color in the
                    array and `cmax` will correspond to the first
                    color.
                showscale
                    Determines whether or not a colorbar is
                    displayed for this trace.

        Returns
        -------
        plotly.graph_objs.layout.Coloraxis
        """
        return self["coloraxis"]

    @coloraxis.setter
    def coloraxis(self, val):
        self["coloraxis"] = val

    # colorscale
    # ----------
    @property
    def colorscale(self):
        """
        The 'colorscale' property is an instance of Colorscale
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Colorscale`
          - A dict of string/value properties that will be passed
            to the Colorscale constructor

            Supported dict properties:

                diverging
                    Sets the default diverging colorscale. Note
                    that `autocolorscale` must be true for this
                    attribute to work.
                sequential
                    Sets the default sequential colorscale for
                    positive values. Note that `autocolorscale`
                    must be true for this attribute to work.
                sequentialminus
                    Sets the default sequential colorscale for
                    negative values. Note that `autocolorscale`
                    must be true for this attribute to work.

        Returns
        -------
        plotly.graph_objs.layout.Colorscale
        """
        return self["colorscale"]

    @colorscale.setter
    def colorscale(self, val):
        self["colorscale"] = val

    # colorway
    # --------
    @property
    def colorway(self):
        """
        Sets the default trace colors.

        The 'colorway' property is a colorlist that may be specified
        as a tuple, list, one-dimensional numpy array, or pandas Series of valid
        color strings

        Returns
        -------
        list
        """
        return self["colorway"]

    @colorway.setter
    def colorway(self, val):
        self["colorway"] = val

    # computed
    # --------
    @property
    def computed(self):
        """
        Placeholder for exporting automargin-impacting values namely
        `margin.t`, `margin.b`, `margin.l` and `margin.r` in "full-
        json" mode.

        The 'computed' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["computed"]

    @computed.setter
    def computed(self, val):
        self["computed"] = val

    # datarevision
    # ------------
    @property
    def datarevision(self):
        """
        If provided, a changed value tells `Plotly.react` that one or
        more data arrays has changed. This way you can modify arrays
        in-place rather than making a complete new copy for an
        incremental change. If NOT provided, `Plotly.react` assumes
        that data arrays are being treated as immutable, thus any data
        array with a different identity from its predecessor contains
        new data.

        The 'datarevision' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["datarevision"]

    @datarevision.setter
    def datarevision(self, val):
        self["datarevision"] = val

    # dragmode
    # --------
    @property
    def dragmode(self):
        """
        Determines the mode of drag interactions. "select" and "lasso"
        apply only to scatter traces with markers or text. "orbit" and
        "turntable" apply only to 3D scenes.

        The 'dragmode' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['zoom', 'pan', 'select', 'lasso', 'drawclosedpath',
                'drawopenpath', 'drawline', 'drawrect', 'drawcircle',
                'orbit', 'turntable', False]

        Returns
        -------
        Any
        """
        return self["dragmode"]

    @dragmode.setter
    def dragmode(self, val):
        self["dragmode"] = val

    # editrevision
    # ------------
    @property
    def editrevision(self):
        """
        Controls persistence of user-driven changes in `editable: true`
        configuration, other than trace names and axis titles. Defaults
        to `layout.uirevision`.

        The 'editrevision' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["editrevision"]

    @editrevision.setter
    def editrevision(self, val):
        self["editrevision"] = val

    # extendfunnelareacolors
    # ----------------------
    @property
    def extendfunnelareacolors(self):
        """
        If `true`, the funnelarea slice colors (whether given by
        `funnelareacolorway` or inherited from `colorway`) will be
        extended to three times its original length by first repeating
        every color 20% lighter then each color 20% darker. This is
        intended to reduce the likelihood of reusing the same color
        when you have many slices, but you can set `false` to disable.
        Colors provided in the trace, using `marker.colors`, are never
        extended.

        The 'extendfunnelareacolors' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["extendfunnelareacolors"]

    @extendfunnelareacolors.setter
    def extendfunnelareacolors(self, val):
        self["extendfunnelareacolors"] = val

    # extendiciclecolors
    # ------------------
    @property
    def extendiciclecolors(self):
        """
        If `true`, the icicle slice colors (whether given by
        `iciclecolorway` or inherited from `colorway`) will be extended
        to three times its original length by first repeating every
        color 20% lighter then each color 20% darker. This is intended
        to reduce the likelihood of reusing the same color when you
        have many slices, but you can set `false` to disable. Colors
        provided in the trace, using `marker.colors`, are never
        extended.

        The 'extendiciclecolors' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["extendiciclecolors"]

    @extendiciclecolors.setter
    def extendiciclecolors(self, val):
        self["extendiciclecolors"] = val

    # extendpiecolors
    # ---------------
    @property
    def extendpiecolors(self):
        """
        If `true`, the pie slice colors (whether given by `piecolorway`
        or inherited from `colorway`) will be extended to three times
        its original length by first repeating every color 20% lighter
        then each color 20% darker. This is intended to reduce the
        likelihood of reusing the same color when you have many slices,
        but you can set `false` to disable. Colors provided in the
        trace, using `marker.colors`, are never extended.

        The 'extendpiecolors' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["extendpiecolors"]

    @extendpiecolors.setter
    def extendpiecolors(self, val):
        self["extendpiecolors"] = val

    # extendsunburstcolors
    # --------------------
    @property
    def extendsunburstcolors(self):
        """
        If `true`, the sunburst slice colors (whether given by
        `sunburstcolorway` or inherited from `colorway`) will be
        extended to three times its original length by first repeating
        every color 20% lighter then each color 20% darker. This is
        intended to reduce the likelihood of reusing the same color
        when you have many slices, but you can set `false` to disable.
        Colors provided in the trace, using `marker.colors`, are never
        extended.

        The 'extendsunburstcolors' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["extendsunburstcolors"]

    @extendsunburstcolors.setter
    def extendsunburstcolors(self, val):
        self["extendsunburstcolors"] = val

    # extendtreemapcolors
    # -------------------
    @property
    def extendtreemapcolors(self):
        """
        If `true`, the treemap slice colors (whether given by
        `treemapcolorway` or inherited from `colorway`) will be
        extended to three times its original length by first repeating
        every color 20% lighter then each color 20% darker. This is
        intended to reduce the likelihood of reusing the same color
        when you have many slices, but you can set `false` to disable.
        Colors provided in the trace, using `marker.colors`, are never
        extended.

        The 'extendtreemapcolors' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["extendtreemapcolors"]

    @extendtreemapcolors.setter
    def extendtreemapcolors(self, val):
        self["extendtreemapcolors"] = val

    # font
    # ----
    @property
    def font(self):
        """
        Sets the global font. Note that fonts used in traces and other
        layout components inherit from the global font.

        The 'font' property is an instance of Font
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Font`
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
        plotly.graph_objs.layout.Font
        """
        return self["font"]

    @font.setter
    def font(self, val):
        self["font"] = val

    # funnelareacolorway
    # ------------------
    @property
    def funnelareacolorway(self):
        """
        Sets the default funnelarea slice colors. Defaults to the main
        `colorway` used for trace colors. If you specify a new list
        here it can still be extended with lighter and darker colors,
        see `extendfunnelareacolors`.

        The 'funnelareacolorway' property is a colorlist that may be specified
        as a tuple, list, one-dimensional numpy array, or pandas Series of valid
        color strings

        Returns
        -------
        list
        """
        return self["funnelareacolorway"]

    @funnelareacolorway.setter
    def funnelareacolorway(self, val):
        self["funnelareacolorway"] = val

    # funnelgap
    # ---------
    @property
    def funnelgap(self):
        """
        Sets the gap (in plot fraction) between bars of adjacent
        location coordinates.

        The 'funnelgap' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["funnelgap"]

    @funnelgap.setter
    def funnelgap(self, val):
        self["funnelgap"] = val

    # funnelgroupgap
    # --------------
    @property
    def funnelgroupgap(self):
        """
        Sets the gap (in plot fraction) between bars of the same
        location coordinate.

        The 'funnelgroupgap' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["funnelgroupgap"]

    @funnelgroupgap.setter
    def funnelgroupgap(self, val):
        self["funnelgroupgap"] = val

    # funnelmode
    # ----------
    @property
    def funnelmode(self):
        """
        Determines how bars at the same location coordinate are
        displayed on the graph. With "stack", the bars are stacked on
        top of one another With "group", the bars are plotted next to
        one another centered around the shared location. With
        "overlay", the bars are plotted over one another, you might
        need to reduce "opacity" to see multiple bars.

        The 'funnelmode' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['stack', 'group', 'overlay']

        Returns
        -------
        Any
        """
        return self["funnelmode"]

    @funnelmode.setter
    def funnelmode(self, val):
        self["funnelmode"] = val

    # geo
    # ---
    @property
    def geo(self):
        """
        The 'geo' property is an instance of Geo
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Geo`
          - A dict of string/value properties that will be passed
            to the Geo constructor

            Supported dict properties:

                bgcolor
                    Set the background color of the map
                center
                    :class:`plotly.graph_objects.layout.geo.Center`
                    instance or dict with compatible properties
                coastlinecolor
                    Sets the coastline color.
                coastlinewidth
                    Sets the coastline stroke width (in px).
                countrycolor
                    Sets line color of the country boundaries.
                countrywidth
                    Sets line width (in px) of the country
                    boundaries.
                domain
                    :class:`plotly.graph_objects.layout.geo.Domain`
                    instance or dict with compatible properties
                fitbounds
                    Determines if this subplot's view settings are
                    auto-computed to fit trace data. On scoped
                    maps, setting `fitbounds` leads to `center.lon`
                    and `center.lat` getting auto-filled. On maps
                    with a non-clipped projection, setting
                    `fitbounds` leads to `center.lon`,
                    `center.lat`, and `projection.rotation.lon`
                    getting auto-filled. On maps with a clipped
                    projection, setting `fitbounds` leads to
                    `center.lon`, `center.lat`,
                    `projection.rotation.lon`,
                    `projection.rotation.lat`, `lonaxis.range` and
                    `lonaxis.range` getting auto-filled. If
                    "locations", only the trace's visible locations
                    are considered in the `fitbounds` computations.
                    If "geojson", the entire trace input `geojson`
                    (if provided) is considered in the `fitbounds`
                    computations, Defaults to False.
                framecolor
                    Sets the color the frame.
                framewidth
                    Sets the stroke width (in px) of the frame.
                lakecolor
                    Sets the color of the lakes.
                landcolor
                    Sets the land mass color.
                lataxis
                    :class:`plotly.graph_objects.layout.geo.Lataxis
                    ` instance or dict with compatible properties
                lonaxis
                    :class:`plotly.graph_objects.layout.geo.Lonaxis
                    ` instance or dict with compatible properties
                oceancolor
                    Sets the ocean color
                projection
                    :class:`plotly.graph_objects.layout.geo.Project
                    ion` instance or dict with compatible
                    properties
                resolution
                    Sets the resolution of the base layers. The
                    values have units of km/mm e.g. 110 corresponds
                    to a scale ratio of 1:110,000,000.
                rivercolor
                    Sets color of the rivers.
                riverwidth
                    Sets the stroke width (in px) of the rivers.
                scope
                    Set the scope of the map.
                showcoastlines
                    Sets whether or not the coastlines are drawn.
                showcountries
                    Sets whether or not country boundaries are
                    drawn.
                showframe
                    Sets whether or not a frame is drawn around the
                    map.
                showlakes
                    Sets whether or not lakes are drawn.
                showland
                    Sets whether or not land masses are filled in
                    color.
                showocean
                    Sets whether or not oceans are filled in color.
                showrivers
                    Sets whether or not rivers are drawn.
                showsubunits
                    Sets whether or not boundaries of subunits
                    within countries (e.g. states, provinces) are
                    drawn.
                subunitcolor
                    Sets the color of the subunits boundaries.
                subunitwidth
                    Sets the stroke width (in px) of the subunits
                    boundaries.
                uirevision
                    Controls persistence of user-driven changes in
                    the view (projection and center). Defaults to
                    `layout.uirevision`.
                visible
                    Sets the default visibility of the base layers.

        Returns
        -------
        plotly.graph_objs.layout.Geo
        """
        return self["geo"]

    @geo.setter
    def geo(self, val):
        self["geo"] = val

    # grid
    # ----
    @property
    def grid(self):
        """
        The 'grid' property is an instance of Grid
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Grid`
          - A dict of string/value properties that will be passed
            to the Grid constructor

            Supported dict properties:

                columns
                    The number of columns in the grid. If you
                    provide a 2D `subplots` array, the length of
                    its longest row is used as the default. If you
                    give an `xaxes` array, its length is used as
                    the default. But it's also possible to have a
                    different length, if you want to leave a row at
                    the end for non-cartesian subplots.
                domain
                    :class:`plotly.graph_objects.layout.grid.Domain
                    ` instance or dict with compatible properties
                pattern
                    If no `subplots`, `xaxes`, or `yaxes` are given
                    but we do have `rows` and `columns`, we can
                    generate defaults using consecutive axis IDs,
                    in two ways: "coupled" gives one x axis per
                    column and one y axis per row. "independent"
                    uses a new xy pair for each cell, left-to-right
                    across each row then iterating rows according
                    to `roworder`.
                roworder
                    Is the first row the top or the bottom? Note
                    that columns are always enumerated from left to
                    right.
                rows
                    The number of rows in the grid. If you provide
                    a 2D `subplots` array or a `yaxes` array, its
                    length is used as the default. But it's also
                    possible to have a different length, if you
                    want to leave a row at the end for non-
                    cartesian subplots.
                subplots
                    Used for freeform grids, where some axes may be
                    shared across subplots but others are not. Each
                    entry should be a cartesian subplot id, like
                    "xy" or "x3y2", or "" to leave that cell empty.
                    You may reuse x axes within the same column,
                    and y axes within the same row. Non-cartesian
                    subplots and traces that support `domain` can
                    place themselves in this grid separately using
                    the `gridcell` attribute.
                xaxes
                    Used with `yaxes` when the x and y axes are
                    shared across columns and rows. Each entry
                    should be an x axis id like "x", "x2", etc., or
                    "" to not put an x axis in that column. Entries
                    other than "" must be unique. Ignored if
                    `subplots` is present. If missing but `yaxes`
                    is present, will generate consecutive IDs.
                xgap
                    Horizontal space between grid cells, expressed
                    as a fraction of the total width available to
                    one cell. Defaults to 0.1 for coupled-axes
                    grids and 0.2 for independent grids.
                xside
                    Sets where the x axis labels and titles go.
                    "bottom" means the very bottom of the grid.
                    "bottom plot" is the lowest plot that each x
                    axis is used in. "top" and "top plot" are
                    similar.
                yaxes
                    Used with `yaxes` when the x and y axes are
                    shared across columns and rows. Each entry
                    should be an y axis id like "y", "y2", etc., or
                    "" to not put a y axis in that row. Entries
                    other than "" must be unique. Ignored if
                    `subplots` is present. If missing but `xaxes`
                    is present, will generate consecutive IDs.
                ygap
                    Vertical space between grid cells, expressed as
                    a fraction of the total height available to one
                    cell. Defaults to 0.1 for coupled-axes grids
                    and 0.3 for independent grids.
                yside
                    Sets where the y axis labels and titles go.
                    "left" means the very left edge of the grid.
                    *left plot* is the leftmost plot that each y
                    axis is used in. "right" and *right plot* are
                    similar.

        Returns
        -------
        plotly.graph_objs.layout.Grid
        """
        return self["grid"]

    @grid.setter
    def grid(self, val):
        self["grid"] = val

    # height
    # ------
    @property
    def height(self):
        """
        Sets the plot's height (in px).

        The 'height' property is a number and may be specified as:
          - An int or float in the interval [10, inf]

        Returns
        -------
        int|float
        """
        return self["height"]

    @height.setter
    def height(self, val):
        self["height"] = val

    # hiddenlabels
    # ------------
    @property
    def hiddenlabels(self):
        """
        hiddenlabels is the funnelarea & pie chart analog of
        visible:'legendonly' but it can contain many labels, and can
        simultaneously hide slices from several pies/funnelarea charts

        The 'hiddenlabels' property is an array that may be specified as a tuple,
        list, numpy array, or pandas Series

        Returns
        -------
        numpy.ndarray
        """
        return self["hiddenlabels"]

    @hiddenlabels.setter
    def hiddenlabels(self, val):
        self["hiddenlabels"] = val

    # hiddenlabelssrc
    # ---------------
    @property
    def hiddenlabelssrc(self):
        """
        Sets the source reference on Chart Studio Cloud for
        `hiddenlabels`.

        The 'hiddenlabelssrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["hiddenlabelssrc"]

    @hiddenlabelssrc.setter
    def hiddenlabelssrc(self, val):
        self["hiddenlabelssrc"] = val

    # hidesources
    # -----------
    @property
    def hidesources(self):
        """
        Determines whether or not a text link citing the data source is
        placed at the bottom-right cored of the figure. Has only an
        effect only on graphs that have been generated via forked
        graphs from the Chart Studio Cloud (at https://chart-
        studio.plotly.com or on-premise).

        The 'hidesources' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["hidesources"]

    @hidesources.setter
    def hidesources(self, val):
        self["hidesources"] = val

    # hoverdistance
    # -------------
    @property
    def hoverdistance(self):
        """
        Sets the default distance (in pixels) to look for data to add
        hover labels (-1 means no cutoff, 0 means no looking for data).
        This is only a real distance for hovering on point-like
        objects, like scatter points. For area-like objects (bars,
        scatter fills, etc) hovering is on inside the area and off
        outside, but these objects will not supersede hover on point-
        like objects in case of conflict.

        The 'hoverdistance' property is a integer and may be specified as:
          - An int (or float that will be cast to an int)
            in the interval [-1, 9223372036854775807]

        Returns
        -------
        int
        """
        return self["hoverdistance"]

    @hoverdistance.setter
    def hoverdistance(self, val):
        self["hoverdistance"] = val

    # hoverlabel
    # ----------
    @property
    def hoverlabel(self):
        """
        The 'hoverlabel' property is an instance of Hoverlabel
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Hoverlabel`
          - A dict of string/value properties that will be passed
            to the Hoverlabel constructor

            Supported dict properties:

                align
                    Sets the horizontal alignment of the text
                    content within hover label box. Has an effect
                    only if the hover label text spans more two or
                    more lines
                bgcolor
                    Sets the background color of all hover labels
                    on graph
                bordercolor
                    Sets the border color of all hover labels on
                    graph.
                font
                    Sets the default hover label font used by all
                    traces on the graph.
                grouptitlefont
                    Sets the font for group titles in hover
                    (unified modes). Defaults to `hoverlabel.font`.
                namelength
                    Sets the default length (in number of
                    characters) of the trace name in the hover
                    labels for all traces. -1 shows the whole name
                    regardless of length. 0-3 shows the first 0-3
                    characters, and an integer >3 will show the
                    whole name if it is less than that many
                    characters, but if it is longer, will truncate
                    to `namelength - 3` characters and add an
                    ellipsis.

        Returns
        -------
        plotly.graph_objs.layout.Hoverlabel
        """
        return self["hoverlabel"]

    @hoverlabel.setter
    def hoverlabel(self, val):
        self["hoverlabel"] = val

    # hovermode
    # ---------
    @property
    def hovermode(self):
        """
        Determines the mode of hover interactions. If "closest", a
        single hoverlabel will appear for the "closest" point within
        the `hoverdistance`. If "x" (or "y"), multiple hoverlabels will
        appear for multiple points at the "closest" x- (or y-)
        coordinate within the `hoverdistance`, with the caveat that no
        more than one hoverlabel will appear per trace. If *x unified*
        (or *y unified*), a single hoverlabel will appear multiple
        points at the closest x- (or y-) coordinate within the
        `hoverdistance` with the caveat that no more than one
        hoverlabel will appear per trace. In this mode, spikelines are
        enabled by default perpendicular to the specified axis. If
        false, hover interactions are disabled.

        The 'hovermode' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['x', 'y', 'closest', False, 'x unified', 'y unified']

        Returns
        -------
        Any
        """
        return self["hovermode"]

    @hovermode.setter
    def hovermode(self, val):
        self["hovermode"] = val

    # hoversubplots
    # -------------
    @property
    def hoversubplots(self):
        """
        Determines expansion of hover effects to other subplots If
        "single" just the axis pair of the primary point is included
        without overlaying subplots. If "overlaying" all subplots using
        the main axis and occupying the same space are included. If
        "axis", also include stacked subplots using the same axis when
        `hovermode` is set to "x", *x unified*, "y" or *y unified*.

        The 'hoversubplots' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['single', 'overlaying', 'axis']

        Returns
        -------
        Any
        """
        return self["hoversubplots"]

    @hoversubplots.setter
    def hoversubplots(self, val):
        self["hoversubplots"] = val

    # iciclecolorway
    # --------------
    @property
    def iciclecolorway(self):
        """
        Sets the default icicle slice colors. Defaults to the main
        `colorway` used for trace colors. If you specify a new list
        here it can still be extended with lighter and darker colors,
        see `extendiciclecolors`.

        The 'iciclecolorway' property is a colorlist that may be specified
        as a tuple, list, one-dimensional numpy array, or pandas Series of valid
        color strings

        Returns
        -------
        list
        """
        return self["iciclecolorway"]

    @iciclecolorway.setter
    def iciclecolorway(self, val):
        self["iciclecolorway"] = val

    # images
    # ------
    @property
    def images(self):
        """
        The 'images' property is a tuple of instances of
        Image that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.Image
          - A list or tuple of dicts of string/value properties that
            will be passed to the Image constructor

            Supported dict properties:

                layer
                    Specifies whether images are drawn below or
                    above traces. When `xref` and `yref` are both
                    set to `paper`, image is drawn below the entire
                    plot area.
                name
                    When used in a template, named items are
                    created in the output figure in addition to any
                    items the figure already has in this array. You
                    can modify these items in the output figure by
                    making your own item with `templateitemname`
                    matching this `name` alongside your
                    modifications (including `visible: false` or
                    `enabled: false` to hide it). Has no effect
                    outside of a template.
                opacity
                    Sets the opacity of the image.
                sizex
                    Sets the image container size horizontally. The
                    image will be sized based on the `position`
                    value. When `xref` is set to `paper`, units are
                    sized relative to the plot width. When `xref`
                    ends with ` domain`, units are sized relative
                    to the axis width.
                sizey
                    Sets the image container size vertically. The
                    image will be sized based on the `position`
                    value. When `yref` is set to `paper`, units are
                    sized relative to the plot height. When `yref`
                    ends with ` domain`, units are sized relative
                    to the axis height.
                sizing
                    Specifies which dimension of the image to
                    constrain.
                source
                    Specifies the URL of the image to be used. The
                    URL must be accessible from the domain where
                    the plot code is run, and can be either
                    relative or absolute.
                templateitemname
                    Used to refer to a named item in this array in
                    the template. Named items from the template
                    will be created even without a matching item in
                    the input figure, but you can modify one by
                    making an item with `templateitemname` matching
                    its `name`, alongside your modifications
                    (including `visible: false` or `enabled: false`
                    to hide it). If there is no template or no
                    matching item, this item will be hidden unless
                    you explicitly show it with `visible: true`.
                visible
                    Determines whether or not this image is
                    visible.
                x
                    Sets the image's x position. When `xref` is set
                    to `paper`, units are sized relative to the
                    plot height. See `xref` for more info
                xanchor
                    Sets the anchor for the x position
                xref
                    Sets the images's x coordinate axis. If set to
                    a x axis id (e.g. "x" or "x2"), the `x`
                    position refers to a x coordinate. If set to
                    "paper", the `x` position refers to the
                    distance from the left of the plotting area in
                    normalized coordinates where 0 (1) corresponds
                    to the left (right). If set to a x axis ID
                    followed by "domain" (separated by a space),
                    the position behaves like for "paper", but
                    refers to the distance in fractions of the
                    domain length from the left of the domain of
                    that axis: e.g., *x2 domain* refers to the
                    domain of the second x  axis and a x position
                    of 0.5 refers to the point between the left and
                    the right of the domain of the second x axis.
                y
                    Sets the image's y position. When `yref` is set
                    to `paper`, units are sized relative to the
                    plot height. See `yref` for more info
                yanchor
                    Sets the anchor for the y position.
                yref
                    Sets the images's y coordinate axis. If set to
                    a y axis id (e.g. "y" or "y2"), the `y`
                    position refers to a y coordinate. If set to
                    "paper", the `y` position refers to the
                    distance from the bottom of the plotting area
                    in normalized coordinates where 0 (1)
                    corresponds to the bottom (top). If set to a y
                    axis ID followed by "domain" (separated by a
                    space), the position behaves like for "paper",
                    but refers to the distance in fractions of the
                    domain length from the bottom of the domain of
                    that axis: e.g., *y2 domain* refers to the
                    domain of the second y  axis and a y position
                    of 0.5 refers to the point between the bottom
                    and the top of the domain of the second y axis.

        Returns
        -------
        tuple[plotly.graph_objs.layout.Image]
        """
        return self["images"]

    @images.setter
    def images(self, val):
        self["images"] = val

    # imagedefaults
    # -------------
    @property
    def imagedefaults(self):
        """
        When used in a template (as
        layout.template.layout.imagedefaults), sets the default
        property values to use for elements of layout.images

        The 'imagedefaults' property is an instance of Image
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Image`
          - A dict of string/value properties that will be passed
            to the Image constructor

            Supported dict properties:

        Returns
        -------
        plotly.graph_objs.layout.Image
        """
        return self["imagedefaults"]

    @imagedefaults.setter
    def imagedefaults(self, val):
        self["imagedefaults"] = val

    # legend
    # ------
    @property
    def legend(self):
        """
        The 'legend' property is an instance of Legend
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Legend`
          - A dict of string/value properties that will be passed
            to the Legend constructor

            Supported dict properties:

                bgcolor
                    Sets the legend background color. Defaults to
                    `layout.paper_bgcolor`.
                bordercolor
                    Sets the color of the border enclosing the
                    legend.
                borderwidth
                    Sets the width (in px) of the border enclosing
                    the legend.
                entrywidth
                    Sets the width (in px or fraction) of the
                    legend. Use 0 to size the entry based on the
                    text width, when `entrywidthmode` is set to
                    "pixels".
                entrywidthmode
                    Determines what entrywidth means.
                font
                    Sets the font used to text the legend items.
                groupclick
                    Determines the behavior on legend group item
                    click. "toggleitem" toggles the visibility of
                    the individual item clicked on the graph.
                    "togglegroup" toggles the visibility of all
                    items in the same legendgroup as the item
                    clicked on the graph.
                grouptitlefont
                    Sets the font for group titles in legend.
                    Defaults to `legend.font` with its size
                    increased about 10%.
                indentation
                    Sets the indentation (in px) of the legend
                    entries.
                itemclick
                    Determines the behavior on legend item click.
                    "toggle" toggles the visibility of the item
                    clicked on the graph. "toggleothers" makes the
                    clicked item the sole visible item on the
                    graph. False disables legend item click
                    interactions.
                itemdoubleclick
                    Determines the behavior on legend item double-
                    click. "toggle" toggles the visibility of the
                    item clicked on the graph. "toggleothers" makes
                    the clicked item the sole visible item on the
                    graph. False disables legend item double-click
                    interactions.
                itemsizing
                    Determines if the legend items symbols scale
                    with their corresponding "trace" attributes or
                    remain "constant" independent of the symbol
                    size on the graph.
                itemwidth
                    Sets the width (in px) of the legend item
                    symbols (the part other than the title.text).
                orientation
                    Sets the orientation of the legend.
                title
                    :class:`plotly.graph_objects.layout.legend.Titl
                    e` instance or dict with compatible properties
                tracegroupgap
                    Sets the amount of vertical space (in px)
                    between legend groups.
                traceorder
                    Determines the order at which the legend items
                    are displayed. If "normal", the items are
                    displayed top-to-bottom in the same order as
                    the input data. If "reversed", the items are
                    displayed in the opposite order as "normal". If
                    "grouped", the items are displayed in groups
                    (when a trace `legendgroup` is provided). if
                    "grouped+reversed", the items are displayed in
                    the opposite order as "grouped".
                uirevision
                    Controls persistence of legend-driven changes
                    in trace and pie label visibility. Defaults to
                    `layout.uirevision`.
                valign
                    Sets the vertical alignment of the symbols with
                    respect to their associated text.
                visible
                    Determines whether or not this legend is
                    visible.
                x
                    Sets the x position with respect to `xref` (in
                    normalized coordinates) of the legend. When
                    `xref` is "paper", defaults to 1.02 for
                    vertical legends and defaults to 0 for
                    horizontal legends. When `xref` is "container",
                    defaults to 1 for vertical legends and defaults
                    to 0 for horizontal legends. Must be between 0
                    and 1 if `xref` is "container". and between
                    "-2" and 3 if `xref` is "paper".
                xanchor
                    Sets the legend's horizontal position anchor.
                    This anchor binds the `x` position to the
                    "left", "center" or "right" of the legend.
                    Value "auto" anchors legends to the right for
                    `x` values greater than or equal to 2/3,
                    anchors legends to the left for `x` values less
                    than or equal to 1/3 and anchors legends with
                    respect to their center otherwise.
                xref
                    Sets the container `x` refers to. "container"
                    spans the entire `width` of the plot. "paper"
                    refers to the width of the plotting area only.
                y
                    Sets the y position with respect to `yref` (in
                    normalized coordinates) of the legend. When
                    `yref` is "paper", defaults to 1 for vertical
                    legends, defaults to "-0.1" for horizontal
                    legends on graphs w/o range sliders and
                    defaults to 1.1 for horizontal legends on graph
                    with one or multiple range sliders. When `yref`
                    is "container", defaults to 1. Must be between
                    0 and 1 if `yref` is "container" and between
                    "-2" and 3 if `yref` is "paper".
                yanchor
                    Sets the legend's vertical position anchor This
                    anchor binds the `y` position to the "top",
                    "middle" or "bottom" of the legend. Value
                    "auto" anchors legends at their bottom for `y`
                    values less than or equal to 1/3, anchors
                    legends to at their top for `y` values greater
                    than or equal to 2/3 and anchors legends with
                    respect to their middle otherwise.
                yref
                    Sets the container `y` refers to. "container"
                    spans the entire `height` of the plot. "paper"
                    refers to the height of the plotting area only.

        Returns
        -------
        plotly.graph_objs.layout.Legend
        """
        return self["legend"]

    @legend.setter
    def legend(self, val):
        self["legend"] = val

    # map
    # ---
    @property
    def map(self):
        """
        The 'map' property is an instance of Map
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Map`
          - A dict of string/value properties that will be passed
            to the Map constructor

            Supported dict properties:

                bearing
                    Sets the bearing angle of the map in degrees
                    counter-clockwise from North (map.bearing).
                bounds
                    :class:`plotly.graph_objects.layout.map.Bounds`
                    instance or dict with compatible properties
                center
                    :class:`plotly.graph_objects.layout.map.Center`
                    instance or dict with compatible properties
                domain
                    :class:`plotly.graph_objects.layout.map.Domain`
                    instance or dict with compatible properties
                layers
                    A tuple of
                    :class:`plotly.graph_objects.layout.map.Layer`
                    instances or dicts with compatible properties
                layerdefaults
                    When used in a template (as
                    layout.template.layout.map.layerdefaults), sets
                    the default property values to use for elements
                    of layout.map.layers
                pitch
                    Sets the pitch angle of the map (in degrees,
                    where 0 means perpendicular to the surface of
                    the map) (map.pitch).
                style
                    Defines the map layers that are rendered by
                    default below the trace layers defined in
                    `data`, which are themselves by default
                    rendered below the layers defined in
                    `layout.map.layers`.  These layers can be
                    defined either explicitly as a Map Style object
                    which can contain multiple layer definitions
                    that load data from any public or private Tile
                    Map Service (TMS or XYZ) or Web Map Service
                    (WMS) or implicitly by using one of the built-
                    in style objects which use WMSes or by using a
                    custom style URL  Map Style objects are of the
                    form described in the MapLibre GL JS
                    documentation available at
                    https://maplibre.org/maplibre-style-spec/  The
                    built-in plotly.js styles objects are: basic,
                    carto-darkmatter, carto-darkmatter-nolabels,
                    carto-positron, carto-positron-nolabels, carto-
                    voyager, carto-voyager-nolabels, dark, light,
                    open-street-map, outdoors, satellite,
                    satellite-streets, streets, white-bg.
                uirevision
                    Controls persistence of user-driven changes in
                    the view: `center`, `zoom`, `bearing`, `pitch`.
                    Defaults to `layout.uirevision`.
                zoom
                    Sets the zoom level of the map (map.zoom).

        Returns
        -------
        plotly.graph_objs.layout.Map
        """
        return self["map"]

    @map.setter
    def map(self, val):
        self["map"] = val

    # mapbox
    # ------
    @property
    def mapbox(self):
        """
        The 'mapbox' property is an instance of Mapbox
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Mapbox`
          - A dict of string/value properties that will be passed
            to the Mapbox constructor

            Supported dict properties:

                accesstoken
                    Sets the mapbox access token to be used for
                    this mapbox map. Alternatively, the mapbox
                    access token can be set in the configuration
                    options under `mapboxAccessToken`. Note that
                    accessToken are only required when `style` (e.g
                    with values : basic, streets, outdoors, light,
                    dark, satellite, satellite-streets ) and/or a
                    layout layer references the Mapbox server.
                bearing
                    Sets the bearing angle of the map in degrees
                    counter-clockwise from North (mapbox.bearing).
                bounds
                    :class:`plotly.graph_objects.layout.mapbox.Boun
                    ds` instance or dict with compatible properties
                center
                    :class:`plotly.graph_objects.layout.mapbox.Cent
                    er` instance or dict with compatible properties
                domain
                    :class:`plotly.graph_objects.layout.mapbox.Doma
                    in` instance or dict with compatible properties
                layers
                    A tuple of :class:`plotly.graph_objects.layout.
                    mapbox.Layer` instances or dicts with
                    compatible properties
                layerdefaults
                    When used in a template (as
                    layout.template.layout.mapbox.layerdefaults),
                    sets the default property values to use for
                    elements of layout.mapbox.layers
                pitch
                    Sets the pitch angle of the map (in degrees,
                    where 0 means perpendicular to the surface of
                    the map) (mapbox.pitch).
                style
                    Defines the map layers that are rendered by
                    default below the trace layers defined in
                    `data`, which are themselves by default
                    rendered below the layers defined in
                    `layout.mapbox.layers`.  These layers can be
                    defined either explicitly as a Mapbox Style
                    object which can contain multiple layer
                    definitions that load data from any public or
                    private Tile Map Service (TMS or XYZ) or Web
                    Map Service (WMS) or implicitly by using one of
                    the built-in style objects which use WMSes
                    which do not require any access tokens, or by
                    using a default Mapbox style or custom Mapbox
                    style URL, both of which require a Mapbox
                    access token  Note that Mapbox access token can
                    be set in the `accesstoken` attribute or in the
                    `mapboxAccessToken` config option.  Mapbox
                    Style objects are of the form described in the
                    Mapbox GL JS documentation available at
                    https://docs.mapbox.com/mapbox-gl-js/style-spec
                    The built-in plotly.js styles objects are:
                    carto-darkmatter, carto-positron, open-street-
                    map, stamen-terrain, stamen-toner, stamen-
                    watercolor, white-bg  The built-in Mapbox
                    styles are: basic, streets, outdoors, light,
                    dark, satellite, satellite-streets  Mapbox
                    style URLs are of the form:
                    mapbox://mapbox.mapbox-<name>-<version>
                uirevision
                    Controls persistence of user-driven changes in
                    the view: `center`, `zoom`, `bearing`, `pitch`.
                    Defaults to `layout.uirevision`.
                zoom
                    Sets the zoom level of the map (mapbox.zoom).

        Returns
        -------
        plotly.graph_objs.layout.Mapbox
        """
        return self["mapbox"]

    @mapbox.setter
    def mapbox(self, val):
        self["mapbox"] = val

    # margin
    # ------
    @property
    def margin(self):
        """
        The 'margin' property is an instance of Margin
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Margin`
          - A dict of string/value properties that will be passed
            to the Margin constructor

            Supported dict properties:

                autoexpand
                    Turns on/off margin expansion computations.
                    Legends, colorbars, updatemenus, sliders, axis
                    rangeselector and rangeslider are allowed to
                    push the margins by defaults.
                b
                    Sets the bottom margin (in px).
                l
                    Sets the left margin (in px).
                pad
                    Sets the amount of padding (in px) between the
                    plotting area and the axis lines
                r
                    Sets the right margin (in px).
                t
                    Sets the top margin (in px).

        Returns
        -------
        plotly.graph_objs.layout.Margin
        """
        return self["margin"]

    @margin.setter
    def margin(self, val):
        self["margin"] = val

    # meta
    # ----
    @property
    def meta(self):
        """
        Assigns extra meta information that can be used in various
        `text` attributes. Attributes such as the graph, axis and
        colorbar `title.text`, annotation `text` `trace.name` in legend
        items, `rangeselector`, `updatemenus` and `sliders` `label`
        text all support `meta`. One can access `meta` fields using
        template strings: `%{meta[i]}` where `i` is the index of the
        `meta` item in question. `meta` can also be an object for
        example `{key: value}` which can be accessed %{meta[key]}.

        The 'meta' property accepts values of any type

        Returns
        -------
        Any|numpy.ndarray
        """
        return self["meta"]

    @meta.setter
    def meta(self, val):
        self["meta"] = val

    # metasrc
    # -------
    @property
    def metasrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `meta`.

        The 'metasrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["metasrc"]

    @metasrc.setter
    def metasrc(self, val):
        self["metasrc"] = val

    # minreducedheight
    # ----------------
    @property
    def minreducedheight(self):
        """
        Minimum height of the plot with margin.automargin applied (in
        px)

        The 'minreducedheight' property is a number and may be specified as:
          - An int or float in the interval [2, inf]

        Returns
        -------
        int|float
        """
        return self["minreducedheight"]

    @minreducedheight.setter
    def minreducedheight(self, val):
        self["minreducedheight"] = val

    # minreducedwidth
    # ---------------
    @property
    def minreducedwidth(self):
        """
        Minimum width of the plot with margin.automargin applied (in
        px)

        The 'minreducedwidth' property is a number and may be specified as:
          - An int or float in the interval [2, inf]

        Returns
        -------
        int|float
        """
        return self["minreducedwidth"]

    @minreducedwidth.setter
    def minreducedwidth(self, val):
        self["minreducedwidth"] = val

    # modebar
    # -------
    @property
    def modebar(self):
        """
        The 'modebar' property is an instance of Modebar
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Modebar`
          - A dict of string/value properties that will be passed
            to the Modebar constructor

            Supported dict properties:

                activecolor
                    Sets the color of the active or hovered on
                    icons in the modebar.
                add
                    Determines which predefined modebar buttons to
                    add. Please note that these buttons will only
                    be shown if they are compatible with all trace
                    types used in a graph. Similar to
                    `config.modeBarButtonsToAdd` option. This may
                    include "v1hovermode", "hoverclosest",
                    "hovercompare", "togglehover",
                    "togglespikelines", "drawline", "drawopenpath",
                    "drawclosedpath", "drawcircle", "drawrect",
                    "eraseshape".
                addsrc
                    Sets the source reference on Chart Studio Cloud
                    for `add`.
                bgcolor
                    Sets the background color of the modebar.
                color
                    Sets the color of the icons in the modebar.
                orientation
                    Sets the orientation of the modebar.
                remove
                    Determines which predefined modebar buttons to
                    remove. Similar to
                    `config.modeBarButtonsToRemove` option. This
                    may include "autoScale2d", "autoscale",
                    "editInChartStudio", "editinchartstudio",
                    "hoverCompareCartesian", "hovercompare",
                    "lasso", "lasso2d", "orbitRotation",
                    "orbitrotation", "pan", "pan2d", "pan3d",
                    "reset", "resetCameraDefault3d",
                    "resetCameraLastSave3d", "resetGeo",
                    "resetSankeyGroup", "resetScale2d",
                    "resetViewMap", "resetViewMapbox",
                    "resetViews", "resetcameradefault",
                    "resetcameralastsave", "resetsankeygroup",
                    "resetscale", "resetview", "resetviews",
                    "select", "select2d", "sendDataToCloud",
                    "senddatatocloud", "tableRotation",
                    "tablerotation", "toImage", "toggleHover",
                    "toggleSpikelines", "togglehover",
                    "togglespikelines", "toimage", "zoom",
                    "zoom2d", "zoom3d", "zoomIn2d", "zoomInGeo",
                    "zoomInMap", "zoomInMapbox", "zoomOut2d",
                    "zoomOutGeo", "zoomOutMap", "zoomOutMapbox",
                    "zoomin", "zoomout".
                removesrc
                    Sets the source reference on Chart Studio Cloud
                    for `remove`.
                uirevision
                    Controls persistence of user-driven changes
                    related to the modebar, including `hovermode`,
                    `dragmode`, and `showspikes` at both the root
                    level and inside subplots. Defaults to
                    `layout.uirevision`.

        Returns
        -------
        plotly.graph_objs.layout.Modebar
        """
        return self["modebar"]

    @modebar.setter
    def modebar(self, val):
        self["modebar"] = val

    # newselection
    # ------------
    @property
    def newselection(self):
        """
        The 'newselection' property is an instance of Newselection
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Newselection`
          - A dict of string/value properties that will be passed
            to the Newselection constructor

            Supported dict properties:

                line
                    :class:`plotly.graph_objects.layout.newselectio
                    n.Line` instance or dict with compatible
                    properties
                mode
                    Describes how a new selection is created. If
                    `immediate`, a new selection is created after
                    first mouse up. If `gradual`, a new selection
                    is not created after first mouse. By adding to
                    and subtracting from the initial selection,
                    this option allows declaring extra outlines of
                    the selection.

        Returns
        -------
        plotly.graph_objs.layout.Newselection
        """
        return self["newselection"]

    @newselection.setter
    def newselection(self, val):
        self["newselection"] = val

    # newshape
    # --------
    @property
    def newshape(self):
        """
        The 'newshape' property is an instance of Newshape
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Newshape`
          - A dict of string/value properties that will be passed
            to the Newshape constructor

            Supported dict properties:

                drawdirection
                    When `dragmode` is set to "drawrect",
                    "drawline" or "drawcircle" this limits the drag
                    to be horizontal, vertical or diagonal. Using
                    "diagonal" there is no limit e.g. in drawing
                    lines in any direction. "ortho" limits the draw
                    to be either horizontal or vertical.
                    "horizontal" allows horizontal extend.
                    "vertical" allows vertical extend.
                fillcolor
                    Sets the color filling new shapes' interior.
                    Please note that if using a fillcolor with
                    alpha greater than half, drag inside the active
                    shape starts moving the shape underneath,
                    otherwise a new shape could be started over.
                fillrule
                    Determines the path's interior. For more info
                    please visit https://developer.mozilla.org/en-
                    US/docs/Web/SVG/Attribute/fill-rule
                label
                    :class:`plotly.graph_objects.layout.newshape.La
                    bel` instance or dict with compatible
                    properties
                layer
                    Specifies whether new shapes are drawn below
                    gridlines ("below"), between gridlines and
                    traces ("between") or above traces ("above").
                legend
                    Sets the reference to a legend to show new
                    shape in. References to these legends are
                    "legend", "legend2", "legend3", etc. Settings
                    for these legends are set in the layout, under
                    `layout.legend`, `layout.legend2`, etc.
                legendgroup
                    Sets the legend group for new shape. Traces and
                    shapes part of the same legend group hide/show
                    at the same time when toggling legend items.
                legendgrouptitle
                    :class:`plotly.graph_objects.layout.newshape.Le
                    gendgrouptitle` instance or dict with
                    compatible properties
                legendrank
                    Sets the legend rank for new shape. Items and
                    groups with smaller ranks are presented on
                    top/left side while with "reversed"
                    `legend.traceorder` they are on bottom/right
                    side. The default legendrank is 1000, so that
                    you can use ranks less than 1000 to place
                    certain items before all unranked items, and
                    ranks greater than 1000 to go after all
                    unranked items.
                legendwidth
                    Sets the width (in px or fraction) of the
                    legend for new shape.
                line
                    :class:`plotly.graph_objects.layout.newshape.Li
                    ne` instance or dict with compatible properties
                name
                    Sets new shape name. The name appears as the
                    legend item.
                opacity
                    Sets the opacity of new shapes.
                showlegend
                    Determines whether or not new shape is shown in
                    the legend.
                visible
                    Determines whether or not new shape is visible.
                    If "legendonly", the shape is not drawn, but
                    can appear as a legend item (provided that the
                    legend itself is visible).

        Returns
        -------
        plotly.graph_objs.layout.Newshape
        """
        return self["newshape"]

    @newshape.setter
    def newshape(self, val):
        self["newshape"] = val

    # paper_bgcolor
    # -------------
    @property
    def paper_bgcolor(self):
        """
        Sets the background color of the paper where the graph is
        drawn.

        The 'paper_bgcolor' property is a color and may be specified as:
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
        return self["paper_bgcolor"]

    @paper_bgcolor.setter
    def paper_bgcolor(self, val):
        self["paper_bgcolor"] = val

    # piecolorway
    # -----------
    @property
    def piecolorway(self):
        """
        Sets the default pie slice colors. Defaults to the main
        `colorway` used for trace colors. If you specify a new list
        here it can still be extended with lighter and darker colors,
        see `extendpiecolors`.

        The 'piecolorway' property is a colorlist that may be specified
        as a tuple, list, one-dimensional numpy array, or pandas Series of valid
        color strings

        Returns
        -------
        list
        """
        return self["piecolorway"]

    @piecolorway.setter
    def piecolorway(self, val):
        self["piecolorway"] = val

    # plot_bgcolor
    # ------------
    @property
    def plot_bgcolor(self):
        """
        Sets the background color of the plotting area in-between x and
        y axes.

        The 'plot_bgcolor' property is a color and may be specified as:
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
        return self["plot_bgcolor"]

    @plot_bgcolor.setter
    def plot_bgcolor(self, val):
        self["plot_bgcolor"] = val

    # polar
    # -----
    @property
    def polar(self):
        """
        The 'polar' property is an instance of Polar
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Polar`
          - A dict of string/value properties that will be passed
            to the Polar constructor

            Supported dict properties:

                angularaxis
                    :class:`plotly.graph_objects.layout.polar.Angul
                    arAxis` instance or dict with compatible
                    properties
                bargap
                    Sets the gap between bars of adjacent location
                    coordinates. Values are unitless, they
                    represent fractions of the minimum difference
                    in bar positions in the data.
                barmode
                    Determines how bars at the same location
                    coordinate are displayed on the graph. With
                    "stack", the bars are stacked on top of one
                    another With "overlay", the bars are plotted
                    over one another, you might need to reduce
                    "opacity" to see multiple bars.
                bgcolor
                    Set the background color of the subplot
                domain
                    :class:`plotly.graph_objects.layout.polar.Domai
                    n` instance or dict with compatible properties
                gridshape
                    Determines if the radial axis grid lines and
                    angular axis line are drawn as "circular"
                    sectors or as "linear" (polygon) sectors. Has
                    an effect only when the angular axis has `type`
                    "category". Note that `radialaxis.angle` is
                    snapped to the angle of the closest vertex when
                    `gridshape` is "circular" (so that radial axis
                    scale is the same as the data scale).
                hole
                    Sets the fraction of the radius to cut out of
                    the polar subplot.
                radialaxis
                    :class:`plotly.graph_objects.layout.polar.Radia
                    lAxis` instance or dict with compatible
                    properties
                sector
                    Sets angular span of this polar subplot with
                    two angles (in degrees). Sector are assumed to
                    be spanned in the counterclockwise direction
                    with 0 corresponding to rightmost limit of the
                    polar subplot.
                uirevision
                    Controls persistence of user-driven changes in
                    axis attributes, if not overridden in the
                    individual axes. Defaults to
                    `layout.uirevision`.

        Returns
        -------
        plotly.graph_objs.layout.Polar
        """
        return self["polar"]

    @polar.setter
    def polar(self, val):
        self["polar"] = val

    # scattergap
    # ----------
    @property
    def scattergap(self):
        """
        Sets the gap (in plot fraction) between scatter points of
        adjacent location coordinates. Defaults to `bargap`.

        The 'scattergap' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["scattergap"]

    @scattergap.setter
    def scattergap(self, val):
        self["scattergap"] = val

    # scattermode
    # -----------
    @property
    def scattermode(self):
        """
        Determines how scatter points at the same location coordinate
        are displayed on the graph. With "group", the scatter points
        are plotted next to one another centered around the shared
        location. With "overlay", the scatter points are plotted over
        one another, you might need to reduce "opacity" to see multiple
        scatter points.

        The 'scattermode' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['group', 'overlay']

        Returns
        -------
        Any
        """
        return self["scattermode"]

    @scattermode.setter
    def scattermode(self, val):
        self["scattermode"] = val

    # scene
    # -----
    @property
    def scene(self):
        """
        The 'scene' property is an instance of Scene
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Scene`
          - A dict of string/value properties that will be passed
            to the Scene constructor

            Supported dict properties:

                annotations
                    A tuple of :class:`plotly.graph_objects.layout.
                    scene.Annotation` instances or dicts with
                    compatible properties
                annotationdefaults
                    When used in a template (as layout.template.lay
                    out.scene.annotationdefaults), sets the default
                    property values to use for elements of
                    layout.scene.annotations
                aspectmode
                    If "cube", this scene's axes are drawn as a
                    cube, regardless of the axes' ranges. If
                    "data", this scene's axes are drawn in
                    proportion with the axes' ranges. If "manual",
                    this scene's axes are drawn in proportion with
                    the input of "aspectratio" (the default
                    behavior if "aspectratio" is provided). If
                    "auto", this scene's axes are drawn using the
                    results of "data" except when one axis is more
                    than four times the size of the two others,
                    where in that case the results of "cube" are
                    used.
                aspectratio
                    Sets this scene's axis aspectratio.
                bgcolor

                camera
                    :class:`plotly.graph_objects.layout.scene.Camer
                    a` instance or dict with compatible properties
                domain
                    :class:`plotly.graph_objects.layout.scene.Domai
                    n` instance or dict with compatible properties
                dragmode
                    Determines the mode of drag interactions for
                    this scene.
                hovermode
                    Determines the mode of hover interactions for
                    this scene.
                uirevision
                    Controls persistence of user-driven changes in
                    camera attributes. Defaults to
                    `layout.uirevision`.
                xaxis
                    :class:`plotly.graph_objects.layout.scene.XAxis
                    ` instance or dict with compatible properties
                yaxis
                    :class:`plotly.graph_objects.layout.scene.YAxis
                    ` instance or dict with compatible properties
                zaxis
                    :class:`plotly.graph_objects.layout.scene.ZAxis
                    ` instance or dict with compatible properties

        Returns
        -------
        plotly.graph_objs.layout.Scene
        """
        return self["scene"]

    @scene.setter
    def scene(self, val):
        self["scene"] = val

    # selectdirection
    # ---------------
    @property
    def selectdirection(self):
        """
        When `dragmode` is set to "select", this limits the selection
        of the drag to horizontal, vertical or diagonal. "h" only
        allows horizontal selection, "v" only vertical, "d" only
        diagonal and "any" sets no limit.

        The 'selectdirection' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['h', 'v', 'd', 'any']

        Returns
        -------
        Any
        """
        return self["selectdirection"]

    @selectdirection.setter
    def selectdirection(self, val):
        self["selectdirection"] = val

    # selectionrevision
    # -----------------
    @property
    def selectionrevision(self):
        """
        Controls persistence of user-driven changes in selected points
        from all traces.

        The 'selectionrevision' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["selectionrevision"]

    @selectionrevision.setter
    def selectionrevision(self, val):
        self["selectionrevision"] = val

    # selections
    # ----------
    @property
    def selections(self):
        """
        The 'selections' property is a tuple of instances of
        Selection that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.Selection
          - A list or tuple of dicts of string/value properties that
            will be passed to the Selection constructor

            Supported dict properties:

                line
                    :class:`plotly.graph_objects.layout.selection.L
                    ine` instance or dict with compatible
                    properties
                name
                    When used in a template, named items are
                    created in the output figure in addition to any
                    items the figure already has in this array. You
                    can modify these items in the output figure by
                    making your own item with `templateitemname`
                    matching this `name` alongside your
                    modifications (including `visible: false` or
                    `enabled: false` to hide it). Has no effect
                    outside of a template.
                opacity
                    Sets the opacity of the selection.
                path
                    For `type` "path" - a valid SVG path similar to
                    `shapes.path` in data coordinates. Allowed
                    segments are: M, L and Z.
                templateitemname
                    Used to refer to a named item in this array in
                    the template. Named items from the template
                    will be created even without a matching item in
                    the input figure, but you can modify one by
                    making an item with `templateitemname` matching
                    its `name`, alongside your modifications
                    (including `visible: false` or `enabled: false`
                    to hide it). If there is no template or no
                    matching item, this item will be hidden unless
                    you explicitly show it with `visible: true`.
                type
                    Specifies the selection type to be drawn. If
                    "rect", a rectangle is drawn linking
                    (`x0`,`y0`), (`x1`,`y0`), (`x1`,`y1`) and
                    (`x0`,`y1`). If "path", draw a custom SVG path
                    using `path`.
                x0
                    Sets the selection's starting x position.
                x1
                    Sets the selection's end x position.
                xref
                    Sets the selection's x coordinate axis. If set
                    to a x axis id (e.g. "x" or "x2"), the `x`
                    position refers to a x coordinate. If set to
                    "paper", the `x` position refers to the
                    distance from the left of the plotting area in
                    normalized coordinates where 0 (1) corresponds
                    to the left (right). If set to a x axis ID
                    followed by "domain" (separated by a space),
                    the position behaves like for "paper", but
                    refers to the distance in fractions of the
                    domain length from the left of the domain of
                    that axis: e.g., *x2 domain* refers to the
                    domain of the second x  axis and a x position
                    of 0.5 refers to the point between the left and
                    the right of the domain of the second x axis.
                y0
                    Sets the selection's starting y position.
                y1
                    Sets the selection's end y position.
                yref
                    Sets the selection's x coordinate axis. If set
                    to a y axis id (e.g. "y" or "y2"), the `y`
                    position refers to a y coordinate. If set to
                    "paper", the `y` position refers to the
                    distance from the bottom of the plotting area
                    in normalized coordinates where 0 (1)
                    corresponds to the bottom (top). If set to a y
                    axis ID followed by "domain" (separated by a
                    space), the position behaves like for "paper",
                    but refers to the distance in fractions of the
                    domain length from the bottom of the domain of
                    that axis: e.g., *y2 domain* refers to the
                    domain of the second y  axis and a y position
                    of 0.5 refers to the point between the bottom
                    and the top of the domain of the second y axis.

        Returns
        -------
        tuple[plotly.graph_objs.layout.Selection]
        """
        return self["selections"]

    @selections.setter
    def selections(self, val):
        self["selections"] = val

    # selectiondefaults
    # -----------------
    @property
    def selectiondefaults(self):
        """
        When used in a template (as
        layout.template.layout.selectiondefaults), sets the default
        property values to use for elements of layout.selections

        The 'selectiondefaults' property is an instance of Selection
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Selection`
          - A dict of string/value properties that will be passed
            to the Selection constructor

            Supported dict properties:

        Returns
        -------
        plotly.graph_objs.layout.Selection
        """
        return self["selectiondefaults"]

    @selectiondefaults.setter
    def selectiondefaults(self, val):
        self["selectiondefaults"] = val

    # separators
    # ----------
    @property
    def separators(self):
        """
        Sets the decimal and thousand separators. For example, *. *
        puts a '.' before decimals and a space between thousands. In
        English locales, dflt is ".," but other locales may alter this
        default.

        The 'separators' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["separators"]

    @separators.setter
    def separators(self, val):
        self["separators"] = val

    # shapes
    # ------
    @property
    def shapes(self):
        """
        The 'shapes' property is a tuple of instances of
        Shape that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.Shape
          - A list or tuple of dicts of string/value properties that
            will be passed to the Shape constructor

            Supported dict properties:

                editable
                    Determines whether the shape could be activated
                    for edit or not. Has no effect when the older
                    editable shapes mode is enabled via
                    `config.editable` or
                    `config.edits.shapePosition`.
                fillcolor
                    Sets the color filling the shape's interior.
                    Only applies to closed shapes.
                fillrule
                    Determines which regions of complex paths
                    constitute the interior. For more info please
                    visit https://developer.mozilla.org/en-
                    US/docs/Web/SVG/Attribute/fill-rule
                label
                    :class:`plotly.graph_objects.layout.shape.Label
                    ` instance or dict with compatible properties
                layer
                    Specifies whether shapes are drawn below
                    gridlines ("below"), between gridlines and
                    traces ("between") or above traces ("above").
                legend
                    Sets the reference to a legend to show this
                    shape in. References to these legends are
                    "legend", "legend2", "legend3", etc. Settings
                    for these legends are set in the layout, under
                    `layout.legend`, `layout.legend2`, etc.
                legendgroup
                    Sets the legend group for this shape. Traces
                    and shapes part of the same legend group
                    hide/show at the same time when toggling legend
                    items.
                legendgrouptitle
                    :class:`plotly.graph_objects.layout.shape.Legen
                    dgrouptitle` instance or dict with compatible
                    properties
                legendrank
                    Sets the legend rank for this shape. Items and
                    groups with smaller ranks are presented on
                    top/left side while with "reversed"
                    `legend.traceorder` they are on bottom/right
                    side. The default legendrank is 1000, so that
                    you can use ranks less than 1000 to place
                    certain items before all unranked items, and
                    ranks greater than 1000 to go after all
                    unranked items. When having unranked or equal
                    rank items shapes would be displayed after
                    traces i.e. according to their order in data
                    and layout.
                legendwidth
                    Sets the width (in px or fraction) of the
                    legend for this shape.
                line
                    :class:`plotly.graph_objects.layout.shape.Line`
                    instance or dict with compatible properties
                name
                    When used in a template, named items are
                    created in the output figure in addition to any
                    items the figure already has in this array. You
                    can modify these items in the output figure by
                    making your own item with `templateitemname`
                    matching this `name` alongside your
                    modifications (including `visible: false` or
                    `enabled: false` to hide it). Has no effect
                    outside of a template.
                opacity
                    Sets the opacity of the shape.
                path
                    For `type` "path" - a valid SVG path with the
                    pixel values replaced by data values in
                    `xsizemode`/`ysizemode` being "scaled" and
                    taken unmodified as pixels relative to
                    `xanchor` and `yanchor` in case of "pixel" size
                    mode. There are a few restrictions / quirks
                    only absolute instructions, not relative. So
                    the allowed segments are: M, L, H, V, Q, C, T,
                    S, and Z arcs (A) are not allowed because
                    radius rx and ry are relative. In the future we
                    could consider supporting relative commands,
                    but we would have to decide on how to handle
                    date and log axes. Note that even as is, Q and
                    C Bezier paths that are smooth on linear axes
                    may not be smooth on log, and vice versa. no
                    chained "polybezier" commands - specify the
                    segment type for each one. On category axes,
                    values are numbers scaled to the serial numbers
                    of categories because using the categories
                    themselves there would be no way to describe
                    fractional positions On data axes: because
                    space and T are both normal components of path
                    strings, we can't use either to separate date
                    from time parts. Therefore we'll use underscore
                    for this purpose: 2015-02-21_13:45:56.789
                showlegend
                    Determines whether or not this shape is shown
                    in the legend.
                templateitemname
                    Used to refer to a named item in this array in
                    the template. Named items from the template
                    will be created even without a matching item in
                    the input figure, but you can modify one by
                    making an item with `templateitemname` matching
                    its `name`, alongside your modifications
                    (including `visible: false` or `enabled: false`
                    to hide it). If there is no template or no
                    matching item, this item will be hidden unless
                    you explicitly show it with `visible: true`.
                type
                    Specifies the shape type to be drawn. If
                    "line", a line is drawn from (`x0`,`y0`) to
                    (`x1`,`y1`) with respect to the axes' sizing
                    mode. If "circle", a circle is drawn from
                    ((`x0`+`x1`)/2, (`y0`+`y1`)/2)) with radius
                    (|(`x0`+`x1`)/2 - `x0`|, |(`y0`+`y1`)/2
                    -`y0`)|) with respect to the axes' sizing mode.
                    If "rect", a rectangle is drawn linking
                    (`x0`,`y0`), (`x1`,`y0`), (`x1`,`y1`),
                    (`x0`,`y1`), (`x0`,`y0`) with respect to the
                    axes' sizing mode. If "path", draw a custom SVG
                    path using `path`. with respect to the axes'
                    sizing mode.
                visible
                    Determines whether or not this shape is
                    visible. If "legendonly", the shape is not
                    drawn, but can appear as a legend item
                    (provided that the legend itself is visible).
                x0
                    Sets the shape's starting x position. See
                    `type` and `xsizemode` for more info.
                x0shift
                    Shifts `x0` away from the center of the
                    category when `xref` is a "category" or
                    "multicategory" axis. -0.5 corresponds to the
                    start of the category and 0.5 corresponds to
                    the end of the category.
                x1
                    Sets the shape's end x position. See `type` and
                    `xsizemode` for more info.
                x1shift
                    Shifts `x1` away from the center of the
                    category when `xref` is a "category" or
                    "multicategory" axis. -0.5 corresponds to the
                    start of the category and 0.5 corresponds to
                    the end of the category.
                xanchor
                    Only relevant in conjunction with `xsizemode`
                    set to "pixel". Specifies the anchor point on
                    the x axis to which `x0`, `x1` and x
                    coordinates within `path` are relative to. E.g.
                    useful to attach a pixel sized shape to a
                    certain data value. No effect when `xsizemode`
                    not set to "pixel".
                xref
                    Sets the shape's x coordinate axis. If set to a
                    x axis id (e.g. "x" or "x2"), the `x` position
                    refers to a x coordinate. If set to "paper",
                    the `x` position refers to the distance from
                    the left of the plotting area in normalized
                    coordinates where 0 (1) corresponds to the left
                    (right). If set to a x axis ID followed by
                    "domain" (separated by a space), the position
                    behaves like for "paper", but refers to the
                    distance in fractions of the domain length from
                    the left of the domain of that axis: e.g., *x2
                    domain* refers to the domain of the second x
                    axis and a x position of 0.5 refers to the
                    point between the left and the right of the
                    domain of the second x axis.
                xsizemode
                    Sets the shapes's sizing mode along the x axis.
                    If set to "scaled", `x0`, `x1` and x
                    coordinates within `path` refer to data values
                    on the x axis or a fraction of the plot area's
                    width (`xref` set to "paper"). If set to
                    "pixel", `xanchor` specifies the x position in
                    terms of data or plot fraction but `x0`, `x1`
                    and x coordinates within `path` are pixels
                    relative to `xanchor`. This way, the shape can
                    have a fixed width while maintaining a position
                    relative to data or plot fraction.
                y0
                    Sets the shape's starting y position. See
                    `type` and `ysizemode` for more info.
                y0shift
                    Shifts `y0` away from the center of the
                    category when `yref` is a "category" or
                    "multicategory" axis. -0.5 corresponds to the
                    start of the category and 0.5 corresponds to
                    the end of the category.
                y1
                    Sets the shape's end y position. See `type` and
                    `ysizemode` for more info.
                y1shift
                    Shifts `y1` away from the center of the
                    category when `yref` is a "category" or
                    "multicategory" axis. -0.5 corresponds to the
                    start of the category and 0.5 corresponds to
                    the end of the category.
                yanchor
                    Only relevant in conjunction with `ysizemode`
                    set to "pixel". Specifies the anchor point on
                    the y axis to which `y0`, `y1` and y
                    coordinates within `path` are relative to. E.g.
                    useful to attach a pixel sized shape to a
                    certain data value. No effect when `ysizemode`
                    not set to "pixel".
                yref
                    Sets the shape's y coordinate axis. If set to a
                    y axis id (e.g. "y" or "y2"), the `y` position
                    refers to a y coordinate. If set to "paper",
                    the `y` position refers to the distance from
                    the bottom of the plotting area in normalized
                    coordinates where 0 (1) corresponds to the
                    bottom (top). If set to a y axis ID followed by
                    "domain" (separated by a space), the position
                    behaves like for "paper", but refers to the
                    distance in fractions of the domain length from
                    the bottom of the domain of that axis: e.g.,
                    *y2 domain* refers to the domain of the second
                    y  axis and a y position of 0.5 refers to the
                    point between the bottom and the top of the
                    domain of the second y axis.
                ysizemode
                    Sets the shapes's sizing mode along the y axis.
                    If set to "scaled", `y0`, `y1` and y
                    coordinates within `path` refer to data values
                    on the y axis or a fraction of the plot area's
                    height (`yref` set to "paper"). If set to
                    "pixel", `yanchor` specifies the y position in
                    terms of data or plot fraction but `y0`, `y1`
                    and y coordinates within `path` are pixels
                    relative to `yanchor`. This way, the shape can
                    have a fixed height while maintaining a
                    position relative to data or plot fraction.

        Returns
        -------
        tuple[plotly.graph_objs.layout.Shape]
        """
        return self["shapes"]

    @shapes.setter
    def shapes(self, val):
        self["shapes"] = val

    # shapedefaults
    # -------------
    @property
    def shapedefaults(self):
        """
        When used in a template (as
        layout.template.layout.shapedefaults), sets the default
        property values to use for elements of layout.shapes

        The 'shapedefaults' property is an instance of Shape
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Shape`
          - A dict of string/value properties that will be passed
            to the Shape constructor

            Supported dict properties:

        Returns
        -------
        plotly.graph_objs.layout.Shape
        """
        return self["shapedefaults"]

    @shapedefaults.setter
    def shapedefaults(self, val):
        self["shapedefaults"] = val

    # showlegend
    # ----------
    @property
    def showlegend(self):
        """
        Determines whether or not a legend is drawn. Default is `true`
        if there is a trace to show and any of these: a) Two or more
        traces would by default be shown in the legend. b) One pie
        trace is shown in the legend. c) One trace is explicitly given
        with `showlegend: true`.

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

    # sliders
    # -------
    @property
    def sliders(self):
        """
        The 'sliders' property is a tuple of instances of
        Slider that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.Slider
          - A list or tuple of dicts of string/value properties that
            will be passed to the Slider constructor

            Supported dict properties:

                active
                    Determines which button (by index starting from
                    0) is considered active.
                activebgcolor
                    Sets the background color of the slider grip
                    while dragging.
                bgcolor
                    Sets the background color of the slider.
                bordercolor
                    Sets the color of the border enclosing the
                    slider.
                borderwidth
                    Sets the width (in px) of the border enclosing
                    the slider.
                currentvalue
                    :class:`plotly.graph_objects.layout.slider.Curr
                    entvalue` instance or dict with compatible
                    properties
                font
                    Sets the font of the slider step labels.
                len
                    Sets the length of the slider This measure
                    excludes the padding of both ends. That is, the
                    slider's length is this length minus the
                    padding on both ends.
                lenmode
                    Determines whether this slider length is set in
                    units of plot "fraction" or in *pixels. Use
                    `len` to set the value.
                minorticklen
                    Sets the length in pixels of minor step tick
                    marks
                name
                    When used in a template, named items are
                    created in the output figure in addition to any
                    items the figure already has in this array. You
                    can modify these items in the output figure by
                    making your own item with `templateitemname`
                    matching this `name` alongside your
                    modifications (including `visible: false` or
                    `enabled: false` to hide it). Has no effect
                    outside of a template.
                pad
                    Set the padding of the slider component along
                    each side.
                steps
                    A tuple of :class:`plotly.graph_objects.layout.
                    slider.Step` instances or dicts with compatible
                    properties
                stepdefaults
                    When used in a template (as
                    layout.template.layout.slider.stepdefaults),
                    sets the default property values to use for
                    elements of layout.slider.steps
                templateitemname
                    Used to refer to a named item in this array in
                    the template. Named items from the template
                    will be created even without a matching item in
                    the input figure, but you can modify one by
                    making an item with `templateitemname` matching
                    its `name`, alongside your modifications
                    (including `visible: false` or `enabled: false`
                    to hide it). If there is no template or no
                    matching item, this item will be hidden unless
                    you explicitly show it with `visible: true`.
                tickcolor
                    Sets the color of the border enclosing the
                    slider.
                ticklen
                    Sets the length in pixels of step tick marks
                tickwidth
                    Sets the tick width (in px).
                transition
                    :class:`plotly.graph_objects.layout.slider.Tran
                    sition` instance or dict with compatible
                    properties
                visible
                    Determines whether or not the slider is
                    visible.
                x
                    Sets the x position (in normalized coordinates)
                    of the slider.
                xanchor
                    Sets the slider's horizontal position anchor.
                    This anchor binds the `x` position to the
                    "left", "center" or "right" of the range
                    selector.
                y
                    Sets the y position (in normalized coordinates)
                    of the slider.
                yanchor
                    Sets the slider's vertical position anchor This
                    anchor binds the `y` position to the "top",
                    "middle" or "bottom" of the range selector.

        Returns
        -------
        tuple[plotly.graph_objs.layout.Slider]
        """
        return self["sliders"]

    @sliders.setter
    def sliders(self, val):
        self["sliders"] = val

    # sliderdefaults
    # --------------
    @property
    def sliderdefaults(self):
        """
        When used in a template (as
        layout.template.layout.sliderdefaults), sets the default
        property values to use for elements of layout.sliders

        The 'sliderdefaults' property is an instance of Slider
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Slider`
          - A dict of string/value properties that will be passed
            to the Slider constructor

            Supported dict properties:

        Returns
        -------
        plotly.graph_objs.layout.Slider
        """
        return self["sliderdefaults"]

    @sliderdefaults.setter
    def sliderdefaults(self, val):
        self["sliderdefaults"] = val

    # smith
    # -----
    @property
    def smith(self):
        """
        The 'smith' property is an instance of Smith
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Smith`
          - A dict of string/value properties that will be passed
            to the Smith constructor

            Supported dict properties:

                bgcolor
                    Set the background color of the subplot
                domain
                    :class:`plotly.graph_objects.layout.smith.Domai
                    n` instance or dict with compatible properties
                imaginaryaxis
                    :class:`plotly.graph_objects.layout.smith.Imagi
                    naryaxis` instance or dict with compatible
                    properties
                realaxis
                    :class:`plotly.graph_objects.layout.smith.Reala
                    xis` instance or dict with compatible
                    properties

        Returns
        -------
        plotly.graph_objs.layout.Smith
        """
        return self["smith"]

    @smith.setter
    def smith(self, val):
        self["smith"] = val

    # spikedistance
    # -------------
    @property
    def spikedistance(self):
        """
        Sets the default distance (in pixels) to look for data to draw
        spikelines to (-1 means no cutoff, 0 means no looking for
        data). As with hoverdistance, distance does not apply to area-
        like objects. In addition, some objects can be hovered on but
        will not generate spikelines, such as scatter fills.

        The 'spikedistance' property is a integer and may be specified as:
          - An int (or float that will be cast to an int)
            in the interval [-1, 9223372036854775807]

        Returns
        -------
        int
        """
        return self["spikedistance"]

    @spikedistance.setter
    def spikedistance(self, val):
        self["spikedistance"] = val

    # sunburstcolorway
    # ----------------
    @property
    def sunburstcolorway(self):
        """
        Sets the default sunburst slice colors. Defaults to the main
        `colorway` used for trace colors. If you specify a new list
        here it can still be extended with lighter and darker colors,
        see `extendsunburstcolors`.

        The 'sunburstcolorway' property is a colorlist that may be specified
        as a tuple, list, one-dimensional numpy array, or pandas Series of valid
        color strings

        Returns
        -------
        list
        """
        return self["sunburstcolorway"]

    @sunburstcolorway.setter
    def sunburstcolorway(self, val):
        self["sunburstcolorway"] = val

    # template
    # --------
    @property
    def template(self):
        """
        Default attributes to be applied to the plot. This should be a
        dict with format: `{'layout': layoutTemplate, 'data':
        {trace_type: [traceTemplate, ...], ...}}` where
        `layoutTemplate` is a dict matching the structure of
        `figure.layout` and `traceTemplate` is a dict matching the
        structure of the trace with type `trace_type` (e.g. 'scatter').
        Alternatively, this may be specified as an instance of
        plotly.graph_objs.layout.Template.  Trace templates are applied
        cyclically to traces of each type. Container arrays (eg
        `annotations`) have special handling: An object ending in
        `defaults` (eg `annotationdefaults`) is applied to each array
        item. But if an item has a `templateitemname` key we look in
        the template array for an item with matching `name` and apply
        that instead. If no matching `name` is found we mark the item
        invisible. Any named template item not referenced is appended
        to the end of the array, so this can be used to add a watermark
        annotation or a logo image, for example. To omit one of these
        items on the plot, make an item with matching
        `templateitemname` and `visible: false`.

        The 'template' property is an instance of Template
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Template`
          - A dict of string/value properties that will be passed
            to the Template constructor

            Supported dict properties:

                data
                    :class:`plotly.graph_objects.layout.template.Da
                    ta` instance or dict with compatible properties
                layout
                    :class:`plotly.graph_objects.Layout` instance
                    or dict with compatible properties

          - The name of a registered template where current registered templates
            are stored in the plotly.io.templates configuration object. The names
            of all registered templates can be retrieved with:
                >>> import plotly.io as pio
                >>> list(pio.templates)  # doctest: +ELLIPSIS
                ['ggplot2', 'seaborn', 'simple_white', 'plotly', 'plotly_white', ...]

          - A string containing multiple registered template names, joined on '+'
            characters (e.g. 'template1+template2'). In this case the resulting
            template is computed by merging together the collection of registered
            templates

        Returns
        -------
        plotly.graph_objs.layout.Template
        """
        return self["template"]

    @template.setter
    def template(self, val):
        self["template"] = val

    # ternary
    # -------
    @property
    def ternary(self):
        """
        The 'ternary' property is an instance of Ternary
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Ternary`
          - A dict of string/value properties that will be passed
            to the Ternary constructor

            Supported dict properties:

                aaxis
                    :class:`plotly.graph_objects.layout.ternary.Aax
                    is` instance or dict with compatible properties
                baxis
                    :class:`plotly.graph_objects.layout.ternary.Bax
                    is` instance or dict with compatible properties
                bgcolor
                    Set the background color of the subplot
                caxis
                    :class:`plotly.graph_objects.layout.ternary.Cax
                    is` instance or dict with compatible properties
                domain
                    :class:`plotly.graph_objects.layout.ternary.Dom
                    ain` instance or dict with compatible
                    properties
                sum
                    The number each triplet should sum to, and the
                    maximum range of each axis
                uirevision
                    Controls persistence of user-driven changes in
                    axis `min` and `title`, if not overridden in
                    the individual axes. Defaults to
                    `layout.uirevision`.

        Returns
        -------
        plotly.graph_objs.layout.Ternary
        """
        return self["ternary"]

    @ternary.setter
    def ternary(self, val):
        self["ternary"] = val

    # title
    # -----
    @property
    def title(self):
        """
        The 'title' property is an instance of Title
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Title`
          - A dict of string/value properties that will be passed
            to the Title constructor

            Supported dict properties:

                automargin
                    Determines whether the title can automatically
                    push the figure margins. If `yref='paper'` then
                    the margin will expand to ensure that the title
                    doesn’t overlap with the edges of the
                    container. If `yref='container'` then the
                    margins will ensure that the title doesn’t
                    overlap with the plot area, tick labels, and
                    axis titles. If `automargin=true` and the
                    margins need to be expanded, then y will be set
                    to a default 1 and yanchor will be set to an
                    appropriate default to ensure that minimal
                    margin space is needed. Note that when
                    `yref='paper'`, only 1 or 0 are allowed y
                    values. Invalid values will be reset to the
                    default 1.
                font
                    Sets the title font. Note that the title's font
                    used to be customized by the now deprecated
                    `titlefont` attribute.
                pad
                    Sets the padding of the title. Each padding
                    value only applies when the corresponding
                    `xanchor`/`yanchor` value is set accordingly.
                    E.g. for left padding to take effect, `xanchor`
                    must be set to "left". The same rule applies if
                    `xanchor`/`yanchor` is determined
                    automatically. Padding is muted if the
                    respective anchor value is "middle*/*center".
                subtitle
                    :class:`plotly.graph_objects.layout.title.Subti
                    tle` instance or dict with compatible
                    properties
                text
                    Sets the plot's title. Note that before the
                    existence of `title.text`, the title's contents
                    used to be defined as the `title` attribute
                    itself. This behavior has been deprecated.
                x
                    Sets the x position with respect to `xref` in
                    normalized coordinates from 0 (left) to 1
                    (right).
                xanchor
                    Sets the title's horizontal alignment with
                    respect to its x position. "left" means that
                    the title starts at x, "right" means that the
                    title ends at x and "center" means that the
                    title's center is at x. "auto" divides `xref`
                    by three and calculates the `xanchor` value
                    automatically based on the value of `x`.
                xref
                    Sets the container `x` refers to. "container"
                    spans the entire `width` of the plot. "paper"
                    refers to the width of the plotting area only.
                y
                    Sets the y position with respect to `yref` in
                    normalized coordinates from 0 (bottom) to 1
                    (top). "auto" places the baseline of the title
                    onto the vertical center of the top margin.
                yanchor
                    Sets the title's vertical alignment with
                    respect to its y position. "top" means that the
                    title's cap line is at y, "bottom" means that
                    the title's baseline is at y and "middle" means
                    that the title's midline is at y. "auto"
                    divides `yref` by three and calculates the
                    `yanchor` value automatically based on the
                    value of `y`.
                yref
                    Sets the container `y` refers to. "container"
                    spans the entire `height` of the plot. "paper"
                    refers to the height of the plotting area only.

        Returns
        -------
        plotly.graph_objs.layout.Title
        """
        return self["title"]

    @title.setter
    def title(self, val):
        self["title"] = val

    # titlefont
    # ---------
    @property
    def titlefont(self):
        """
        Deprecated: Please use layout.title.font instead. Sets the
        title font. Note that the title's font used to be customized by
        the now deprecated `titlefont` attribute.

        The 'font' property is an instance of Font
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.title.Font`
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

        """
        return self["titlefont"]

    @titlefont.setter
    def titlefont(self, val):
        self["titlefont"] = val

    # transition
    # ----------
    @property
    def transition(self):
        """
        Sets transition options used during Plotly.react updates.

        The 'transition' property is an instance of Transition
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Transition`
          - A dict of string/value properties that will be passed
            to the Transition constructor

            Supported dict properties:

                duration
                    The duration of the transition, in
                    milliseconds. If equal to zero, updates are
                    synchronous.
                easing
                    The easing function used for the transition
                ordering
                    Determines whether the figure's layout or
                    traces smoothly transitions during updates that
                    make both traces and layout change.

        Returns
        -------
        plotly.graph_objs.layout.Transition
        """
        return self["transition"]

    @transition.setter
    def transition(self, val):
        self["transition"] = val

    # treemapcolorway
    # ---------------
    @property
    def treemapcolorway(self):
        """
        Sets the default treemap slice colors. Defaults to the main
        `colorway` used for trace colors. If you specify a new list
        here it can still be extended with lighter and darker colors,
        see `extendtreemapcolors`.

        The 'treemapcolorway' property is a colorlist that may be specified
        as a tuple, list, one-dimensional numpy array, or pandas Series of valid
        color strings

        Returns
        -------
        list
        """
        return self["treemapcolorway"]

    @treemapcolorway.setter
    def treemapcolorway(self, val):
        self["treemapcolorway"] = val

    # uirevision
    # ----------
    @property
    def uirevision(self):
        """
        Used to allow user interactions with the plot to persist after
        `Plotly.react` calls that are unaware of these interactions. If
        `uirevision` is omitted, or if it is given and it changed from
        the previous `Plotly.react` call, the exact new figure is used.
        If `uirevision` is truthy and did NOT change, any attribute
        that has been affected by user interactions and did not receive
        a different value in the new figure will keep the interaction
        value. `layout.uirevision` attribute serves as the default for
        `uirevision` attributes in various sub-containers. For finer
        control you can set these sub-attributes directly. For example,
        if your app separately controls the data on the x and y axes
        you might set `xaxis.uirevision=*time*` and
        `yaxis.uirevision=*cost*`. Then if only the y data is changed,
        you can update `yaxis.uirevision=*quantity*` and the y axis
        range will reset but the x axis range will retain any user-
        driven zoom.

        The 'uirevision' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["uirevision"]

    @uirevision.setter
    def uirevision(self, val):
        self["uirevision"] = val

    # uniformtext
    # -----------
    @property
    def uniformtext(self):
        """
        The 'uniformtext' property is an instance of Uniformtext
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Uniformtext`
          - A dict of string/value properties that will be passed
            to the Uniformtext constructor

            Supported dict properties:

                minsize
                    Sets the minimum text size between traces of
                    the same type.
                mode
                    Determines how the font size for various text
                    elements are uniformed between each trace type.
                    If the computed text sizes were smaller than
                    the minimum size defined by
                    `uniformtext.minsize` using "hide" option hides
                    the text; and using "show" option shows the
                    text without further downscaling. Please note
                    that if the size defined by `minsize` is
                    greater than the font size defined by trace,
                    then the `minsize` is used.

        Returns
        -------
        plotly.graph_objs.layout.Uniformtext
        """
        return self["uniformtext"]

    @uniformtext.setter
    def uniformtext(self, val):
        self["uniformtext"] = val

    # updatemenus
    # -----------
    @property
    def updatemenus(self):
        """
        The 'updatemenus' property is a tuple of instances of
        Updatemenu that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.Updatemenu
          - A list or tuple of dicts of string/value properties that
            will be passed to the Updatemenu constructor

            Supported dict properties:

                active
                    Determines which button (by index starting from
                    0) is considered active.
                bgcolor
                    Sets the background color of the update menu
                    buttons.
                bordercolor
                    Sets the color of the border enclosing the
                    update menu.
                borderwidth
                    Sets the width (in px) of the border enclosing
                    the update menu.
                buttons
                    A tuple of :class:`plotly.graph_objects.layout.
                    updatemenu.Button` instances or dicts with
                    compatible properties
                buttondefaults
                    When used in a template (as layout.template.lay
                    out.updatemenu.buttondefaults), sets the
                    default property values to use for elements of
                    layout.updatemenu.buttons
                direction
                    Determines the direction in which the buttons
                    are laid out, whether in a dropdown menu or a
                    row/column of buttons. For `left` and `up`, the
                    buttons will still appear in left-to-right or
                    top-to-bottom order respectively.
                font
                    Sets the font of the update menu button text.
                name
                    When used in a template, named items are
                    created in the output figure in addition to any
                    items the figure already has in this array. You
                    can modify these items in the output figure by
                    making your own item with `templateitemname`
                    matching this `name` alongside your
                    modifications (including `visible: false` or
                    `enabled: false` to hide it). Has no effect
                    outside of a template.
                pad
                    Sets the padding around the buttons or dropdown
                    menu.
                showactive
                    Highlights active dropdown item or active
                    button if true.
                templateitemname
                    Used to refer to a named item in this array in
                    the template. Named items from the template
                    will be created even without a matching item in
                    the input figure, but you can modify one by
                    making an item with `templateitemname` matching
                    its `name`, alongside your modifications
                    (including `visible: false` or `enabled: false`
                    to hide it). If there is no template or no
                    matching item, this item will be hidden unless
                    you explicitly show it with `visible: true`.
                type
                    Determines whether the buttons are accessible
                    via a dropdown menu or whether the buttons are
                    stacked horizontally or vertically
                visible
                    Determines whether or not the update menu is
                    visible.
                x
                    Sets the x position (in normalized coordinates)
                    of the update menu.
                xanchor
                    Sets the update menu's horizontal position
                    anchor. This anchor binds the `x` position to
                    the "left", "center" or "right" of the range
                    selector.
                y
                    Sets the y position (in normalized coordinates)
                    of the update menu.
                yanchor
                    Sets the update menu's vertical position anchor
                    This anchor binds the `y` position to the
                    "top", "middle" or "bottom" of the range
                    selector.

        Returns
        -------
        tuple[plotly.graph_objs.layout.Updatemenu]
        """
        return self["updatemenus"]

    @updatemenus.setter
    def updatemenus(self, val):
        self["updatemenus"] = val

    # updatemenudefaults
    # ------------------
    @property
    def updatemenudefaults(self):
        """
        When used in a template (as
        layout.template.layout.updatemenudefaults), sets the default
        property values to use for elements of layout.updatemenus

        The 'updatemenudefaults' property is an instance of Updatemenu
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Updatemenu`
          - A dict of string/value properties that will be passed
            to the Updatemenu constructor

            Supported dict properties:

        Returns
        -------
        plotly.graph_objs.layout.Updatemenu
        """
        return self["updatemenudefaults"]

    @updatemenudefaults.setter
    def updatemenudefaults(self, val):
        self["updatemenudefaults"] = val

    # violingap
    # ---------
    @property
    def violingap(self):
        """
        Sets the gap (in plot fraction) between violins of adjacent
        location coordinates. Has no effect on traces that have "width"
        set.

        The 'violingap' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["violingap"]

    @violingap.setter
    def violingap(self, val):
        self["violingap"] = val

    # violingroupgap
    # --------------
    @property
    def violingroupgap(self):
        """
        Sets the gap (in plot fraction) between violins of the same
        location coordinate. Has no effect on traces that have "width"
        set.

        The 'violingroupgap' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["violingroupgap"]

    @violingroupgap.setter
    def violingroupgap(self, val):
        self["violingroupgap"] = val

    # violinmode
    # ----------
    @property
    def violinmode(self):
        """
        Determines how violins at the same location coordinate are
        displayed on the graph. If "group", the violins are plotted
        next to one another centered around the shared location. If
        "overlay", the violins are plotted over one another, you might
        need to set "opacity" to see them multiple violins. Has no
        effect on traces that have "width" set.

        The 'violinmode' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['group', 'overlay']

        Returns
        -------
        Any
        """
        return self["violinmode"]

    @violinmode.setter
    def violinmode(self, val):
        self["violinmode"] = val

    # waterfallgap
    # ------------
    @property
    def waterfallgap(self):
        """
        Sets the gap (in plot fraction) between bars of adjacent
        location coordinates.

        The 'waterfallgap' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["waterfallgap"]

    @waterfallgap.setter
    def waterfallgap(self, val):
        self["waterfallgap"] = val

    # waterfallgroupgap
    # -----------------
    @property
    def waterfallgroupgap(self):
        """
        Sets the gap (in plot fraction) between bars of the same
        location coordinate.

        The 'waterfallgroupgap' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["waterfallgroupgap"]

    @waterfallgroupgap.setter
    def waterfallgroupgap(self, val):
        self["waterfallgroupgap"] = val

    # waterfallmode
    # -------------
    @property
    def waterfallmode(self):
        """
        Determines how bars at the same location coordinate are
        displayed on the graph. With "group", the bars are plotted next
        to one another centered around the shared location. With
        "overlay", the bars are plotted over one another, you might
        need to reduce "opacity" to see multiple bars.

        The 'waterfallmode' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['group', 'overlay']

        Returns
        -------
        Any
        """
        return self["waterfallmode"]

    @waterfallmode.setter
    def waterfallmode(self, val):
        self["waterfallmode"] = val

    # width
    # -----
    @property
    def width(self):
        """
        Sets the plot's width (in px).

        The 'width' property is a number and may be specified as:
          - An int or float in the interval [10, inf]

        Returns
        -------
        int|float
        """
        return self["width"]

    @width.setter
    def width(self, val):
        self["width"] = val

    # xaxis
    # -----
    @property
    def xaxis(self):
        """
        The 'xaxis' property is an instance of XAxis
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.XAxis`
          - A dict of string/value properties that will be passed
            to the XAxis constructor

            Supported dict properties:

                anchor
                    If set to an opposite-letter axis id (e.g.
                    `x2`, `y`), this axis is bound to the
                    corresponding opposite-letter axis. If set to
                    "free", this axis' position is determined by
                    `position`.
                automargin
                    Determines whether long tick labels
                    automatically grow the figure margins.
                autorange
                    Determines whether or not the range of this
                    axis is computed in relation to the input data.
                    See `rangemode` for more info. If `range` is
                    provided and it has a value for both the lower
                    and upper bound, `autorange` is set to False.
                    Using "min" applies autorange only to set the
                    minimum. Using "max" applies autorange only to
                    set the maximum. Using *min reversed* applies
                    autorange only to set the minimum on a reversed
                    axis. Using *max reversed* applies autorange
                    only to set the maximum on a reversed axis.
                    Using "reversed" applies autorange on both ends
                    and reverses the axis direction.
                autorangeoptions
                    :class:`plotly.graph_objects.layout.xaxis.Autor
                    angeoptions` instance or dict with compatible
                    properties
                autotickangles
                    When `tickangle` is set to "auto", it will be
                    set to the first angle in this array that is
                    large enough to prevent label overlap.
                autotypenumbers
                    Using "strict" a numeric string in trace data
                    is not converted to a number. Using *convert
                    types* a numeric string in trace data may be
                    treated as a number during automatic axis
                    `type` detection. Defaults to
                    layout.autotypenumbers.
                calendar
                    Sets the calendar system to use for `range` and
                    `tick0` if this is a date axis. This does not
                    set the calendar for interpreting data on this
                    axis, that's specified in the trace or via the
                    global `layout.calendar`
                categoryarray
                    Sets the order in which categories on this axis
                    appear. Only has an effect if `categoryorder`
                    is set to "array". Used with `categoryorder`.
                categoryarraysrc
                    Sets the source reference on Chart Studio Cloud
                    for `categoryarray`.
                categoryorder
                    Specifies the ordering logic for the case of
                    categorical variables. By default, plotly uses
                    "trace", which specifies the order that is
                    present in the data supplied. Set
                    `categoryorder` to *category ascending* or
                    *category descending* if order should be
                    determined by the alphanumerical order of the
                    category names. Set `categoryorder` to "array"
                    to derive the ordering from the attribute
                    `categoryarray`. If a category is not found in
                    the `categoryarray` array, the sorting behavior
                    for that attribute will be identical to the
                    "trace" mode. The unspecified categories will
                    follow the categories in `categoryarray`. Set
                    `categoryorder` to *total ascending* or *total
                    descending* if order should be determined by
                    the numerical order of the values. Similarly,
                    the order can be determined by the min, max,
                    sum, mean, geometric mean or median of all the
                    values.
                color
                    Sets default for all colors associated with
                    this axis all at once: line, font, tick, and
                    grid colors. Grid color is lightened by
                    blending this with the plot background
                    Individual pieces can override this.
                constrain
                    If this axis needs to be compressed (either due
                    to its own `scaleanchor` and `scaleratio` or
                    those of the other axis), determines how that
                    happens: by increasing the "range", or by
                    decreasing the "domain". Default is "domain"
                    for axes containing image traces, "range"
                    otherwise.
                constraintoward
                    If this axis needs to be compressed (either due
                    to its own `scaleanchor` and `scaleratio` or
                    those of the other axis), determines which
                    direction we push the originally specified plot
                    area. Options are "left", "center" (default),
                    and "right" for x axes, and "top", "middle"
                    (default), and "bottom" for y axes.
                dividercolor
                    Sets the color of the dividers Only has an
                    effect on "multicategory" axes.
                dividerwidth
                    Sets the width (in px) of the dividers Only has
                    an effect on "multicategory" axes.
                domain
                    Sets the domain of this axis (in plot
                    fraction).
                dtick
                    Sets the step in-between ticks on this axis.
                    Use with `tick0`. Must be a positive number, or
                    special strings available to "log" and "date"
                    axes. If the axis `type` is "log", then ticks
                    are set every 10^(n*dtick) where n is the tick
                    number. For example, to set a tick mark at 1,
                    10, 100, 1000, ... set dtick to 1. To set tick
                    marks at 1, 100, 10000, ... set dtick to 2. To
                    set tick marks at 1, 5, 25, 125, 625, 3125, ...
                    set dtick to log_10(5), or 0.69897000433. "log"
                    has several special values; "L<f>", where `f`
                    is a positive number, gives ticks linearly
                    spaced in value (but not position). For example
                    `tick0` = 0.1, `dtick` = "L0.5" will put ticks
                    at 0.1, 0.6, 1.1, 1.6 etc. To show powers of 10
                    plus small digits between, use "D1" (all
                    digits) or "D2" (only 2 and 5). `tick0` is
                    ignored for "D1" and "D2". If the axis `type`
                    is "date", then you must convert the time to
                    milliseconds. For example, to set the interval
                    between ticks to one day, set `dtick` to
                    86400000.0. "date" also has special values
                    "M<n>" gives ticks spaced by a number of
                    months. `n` must be a positive integer. To set
                    ticks on the 15th of every third month, set
                    `tick0` to "2000-01-15" and `dtick` to "M3". To
                    set ticks every 4 years, set `dtick` to "M48"
                exponentformat
                    Determines a formatting rule for the tick
                    exponents. For example, consider the number
                    1,000,000,000. If "none", it appears as
                    1,000,000,000. If "e", 1e+9. If "E", 1E+9. If
                    "power", 1x10^9 (with 9 in a super script). If
                    "SI", 1G. If "B", 1B.
                fixedrange
                    Determines whether or not this axis is zoom-
                    able. If true, then zoom is disabled.
                gridcolor
                    Sets the color of the grid lines.
                griddash
                    Sets the dash style of lines. Set to a dash
                    type string ("solid", "dot", "dash",
                    "longdash", "dashdot", or "longdashdot") or a
                    dash length list in px (eg "5px,10px,2px,2px").
                gridwidth
                    Sets the width (in px) of the grid lines.
                hoverformat
                    Sets the hover text formatting rule using d3
                    formatting mini-languages which are very
                    similar to those in Python. For numbers, see: h
                    ttps://github.com/d3/d3-format/tree/v1.4.5#d3-
                    format. And for dates see:
                    https://github.com/d3/d3-time-
                    format/tree/v2.2.3#locale_format. We add two
                    items to d3's date formatter: "%h" for half of
                    the year as a decimal number as well as "%{n}f"
                    for fractional seconds with n digits. For
                    example, *2016-10-13 09:15:23.456* with
                    tickformat "%H~%M~%S.%2f" would display
                    "09~15~23.46"
                insiderange
                    Could be used to set the desired inside range
                    of this axis (excluding the labels) when
                    `ticklabelposition` of the anchored axis has
                    "inside". Not implemented for axes with `type`
                    "log". This would be ignored when `range` is
                    provided.
                labelalias
                    Replacement text for specific tick or hover
                    labels. For example using {US: 'USA', CA:
                    'Canada'} changes US to USA and CA to Canada.
                    The labels we would have shown must match the
                    keys exactly, after adding any tickprefix or
                    ticksuffix. For negative numbers the minus sign
                    symbol used (U+2212) is wider than the regular
                    ascii dash. That means you need to use −1
                    instead of -1. labelalias can be used with any
                    axis type, and both keys (if needed) and values
                    (if desired) can include html-like tags or
                    MathJax.
                layer
                    Sets the layer on which this axis is displayed.
                    If *above traces*, this axis is displayed above
                    all the subplot's traces If *below traces*,
                    this axis is displayed below all the subplot's
                    traces, but above the grid lines. Useful when
                    used together with scatter-like traces with
                    `cliponaxis` set to False to show markers
                    and/or text nodes above this axis.
                linecolor
                    Sets the axis line color.
                linewidth
                    Sets the width (in px) of the axis line.
                matches
                    If set to another axis id (e.g. `x2`, `y`), the
                    range of this axis will match the range of the
                    corresponding axis in data-coordinates space.
                    Moreover, matching axes share auto-range
                    values, category lists and histogram auto-bins.
                    Note that setting axes simultaneously in both a
                    `scaleanchor` and a `matches` constraint is
                    currently forbidden. Moreover, note that
                    matching axes must have the same `type`.
                maxallowed
                    Determines the maximum range of this axis.
                minallowed
                    Determines the minimum range of this axis.
                minexponent
                    Hide SI prefix for 10^n if |n| is below this
                    number. This only has an effect when
                    `tickformat` is "SI" or "B".
                minor
                    :class:`plotly.graph_objects.layout.xaxis.Minor
                    ` instance or dict with compatible properties
                mirror
                    Determines if the axis lines or/and ticks are
                    mirrored to the opposite side of the plotting
                    area. If True, the axis lines are mirrored. If
                    "ticks", the axis lines and ticks are mirrored.
                    If False, mirroring is disable. If "all", axis
                    lines are mirrored on all shared-axes subplots.
                    If "allticks", axis lines and ticks are
                    mirrored on all shared-axes subplots.
                nticks
                    Specifies the maximum number of ticks for the
                    particular axis. The actual number of ticks
                    will be chosen automatically to be less than or
                    equal to `nticks`. Has an effect only if
                    `tickmode` is set to "auto".
                overlaying
                    If set a same-letter axis id, this axis is
                    overlaid on top of the corresponding same-
                    letter axis, with traces and axes visible for
                    both axes. If False, this axis does not overlay
                    any same-letter axes. In this case, for axes
                    with overlapping domains only the highest-
                    numbered axis will be visible.
                position
                    Sets the position of this axis in the plotting
                    space (in normalized coordinates). Only has an
                    effect if `anchor` is set to "free".
                range
                    Sets the range of this axis. If the axis `type`
                    is "log", then you must take the log of your
                    desired range (e.g. to set the range from 1 to
                    100, set the range from 0 to 2). If the axis
                    `type` is "date", it should be date strings,
                    like date data, though Date objects and unix
                    milliseconds will be accepted and converted to
                    strings. If the axis `type` is "category", it
                    should be numbers, using the scale where each
                    category is assigned a serial number from zero
                    in the order it appears. Leaving either or both
                    elements `null` impacts the default
                    `autorange`.
                rangebreaks
                    A tuple of :class:`plotly.graph_objects.layout.
                    xaxis.Rangebreak` instances or dicts with
                    compatible properties
                rangebreakdefaults
                    When used in a template (as layout.template.lay
                    out.xaxis.rangebreakdefaults), sets the default
                    property values to use for elements of
                    layout.xaxis.rangebreaks
                rangemode
                    If "normal", the range is computed in relation
                    to the extrema of the input data. If *tozero*`,
                    the range extends to 0, regardless of the input
                    data If "nonnegative", the range is non-
                    negative, regardless of the input data. Applies
                    only to linear axes.
                rangeselector
                    :class:`plotly.graph_objects.layout.xaxis.Range
                    selector` instance or dict with compatible
                    properties
                rangeslider
                    :class:`plotly.graph_objects.layout.xaxis.Range
                    slider` instance or dict with compatible
                    properties
                scaleanchor
                    If set to another axis id (e.g. `x2`, `y`), the
                    range of this axis changes together with the
                    range of the corresponding axis such that the
                    scale of pixels per unit is in a constant
                    ratio. Both axes are still zoomable, but when
                    you zoom one, the other will zoom the same
                    amount, keeping a fixed midpoint. `constrain`
                    and `constraintoward` determine how we enforce
                    the constraint. You can chain these, ie `yaxis:
                    {scaleanchor: *x*}, xaxis2: {scaleanchor: *y*}`
                    but you can only link axes of the same `type`.
                    The linked axis can have the opposite letter
                    (to constrain the aspect ratio) or the same
                    letter (to match scales across subplots). Loops
                    (`yaxis: {scaleanchor: *x*}, xaxis:
                    {scaleanchor: *y*}` or longer) are redundant
                    and the last constraint encountered will be
                    ignored to avoid possible inconsistent
                    constraints via `scaleratio`. Note that setting
                    axes simultaneously in both a `scaleanchor` and
                    a `matches` constraint is currently forbidden.
                    Setting `false` allows to remove a default
                    constraint (occasionally, you may need to
                    prevent a default `scaleanchor` constraint from
                    being applied, eg. when having an image trace
                    `yaxis: {scaleanchor: "x"}` is set
                    automatically in order for pixels to be
                    rendered as squares, setting `yaxis:
                    {scaleanchor: false}` allows to remove the
                    constraint).
                scaleratio
                    If this axis is linked to another by
                    `scaleanchor`, this determines the pixel to
                    unit scale ratio. For example, if this value is
                    10, then every unit on this axis spans 10 times
                    the number of pixels as a unit on the linked
                    axis. Use this for example to create an
                    elevation profile where the vertical scale is
                    exaggerated a fixed amount with respect to the
                    horizontal.
                separatethousands
                    If "true", even 4-digit integers are separated
                showdividers
                    Determines whether or not a dividers are drawn
                    between the category levels of this axis. Only
                    has an effect on "multicategory" axes.
                showexponent
                    If "all", all exponents are shown besides their
                    significands. If "first", only the exponent of
                    the first tick is shown. If "last", only the
                    exponent of the last tick is shown. If "none",
                    no exponents appear.
                showgrid
                    Determines whether or not grid lines are drawn.
                    If True, the grid lines are drawn at every tick
                    mark.
                showline
                    Determines whether or not a line bounding this
                    axis is drawn.
                showspikes
                    Determines whether or not spikes (aka
                    droplines) are drawn for this axis. Note: This
                    only takes affect when hovermode = closest
                showticklabels
                    Determines whether or not the tick labels are
                    drawn.
                showtickprefix
                    If "all", all tick labels are displayed with a
                    prefix. If "first", only the first tick is
                    displayed with a prefix. If "last", only the
                    last tick is displayed with a suffix. If
                    "none", tick prefixes are hidden.
                showticksuffix
                    Same as `showtickprefix` but for tick suffixes.
                side
                    Determines whether a x (y) axis is positioned
                    at the "bottom" ("left") or "top" ("right") of
                    the plotting area.
                spikecolor
                    Sets the spike color. If undefined, will use
                    the series color
                spikedash
                    Sets the dash style of lines. Set to a dash
                    type string ("solid", "dot", "dash",
                    "longdash", "dashdot", or "longdashdot") or a
                    dash length list in px (eg "5px,10px,2px,2px").
                spikemode
                    Determines the drawing mode for the spike line
                    If "toaxis", the line is drawn from the data
                    point to the axis the  series is plotted on. If
                    "across", the line is drawn across the entire
                    plot area, and supercedes "toaxis". If
                    "marker", then a marker dot is drawn on the
                    axis the series is plotted on
                spikesnap
                    Determines whether spikelines are stuck to the
                    cursor or to the closest datapoints.
                spikethickness
                    Sets the width (in px) of the zero line.
                tick0
                    Sets the placement of the first tick on this
                    axis. Use with `dtick`. If the axis `type` is
                    "log", then you must take the log of your
                    starting tick (e.g. to set the starting tick to
                    100, set the `tick0` to 2) except when
                    `dtick`=*L<f>* (see `dtick` for more info). If
                    the axis `type` is "date", it should be a date
                    string, like date data. If the axis `type` is
                    "category", it should be a number, using the
                    scale where each category is assigned a serial
                    number from zero in the order it appears.
                tickangle
                    Sets the angle of the tick labels with respect
                    to the horizontal. For example, a `tickangle`
                    of -90 draws the tick labels vertically.
                tickcolor
                    Sets the tick color.
                tickfont
                    Sets the tick font.
                tickformat
                    Sets the tick label formatting rule using d3
                    formatting mini-languages which are very
                    similar to those in Python. For numbers, see: h
                    ttps://github.com/d3/d3-format/tree/v1.4.5#d3-
                    format. And for dates see:
                    https://github.com/d3/d3-time-
                    format/tree/v2.2.3#locale_format. We add two
                    items to d3's date formatter: "%h" for half of
                    the year as a decimal number as well as "%{n}f"
                    for fractional seconds with n digits. For
                    example, *2016-10-13 09:15:23.456* with
                    tickformat "%H~%M~%S.%2f" would display
                    "09~15~23.46"
                tickformatstops
                    A tuple of :class:`plotly.graph_objects.layout.
                    xaxis.Tickformatstop` instances or dicts with
                    compatible properties
                tickformatstopdefaults
                    When used in a template (as layout.template.lay
                    out.xaxis.tickformatstopdefaults), sets the
                    default property values to use for elements of
                    layout.xaxis.tickformatstops
                ticklabelindex
                    Only for axes with `type` "date" or "linear".
                    Instead of drawing the major tick label, draw
                    the label for the minor tick that is n
                    positions away from the major tick. E.g. to
                    always draw the label for the minor tick before
                    each major tick, choose `ticklabelindex` -1.
                    This is useful for date axes with
                    `ticklabelmode` "period" if you want to label
                    the period that ends with each major tick
                    instead of the period that begins there.
                ticklabelindexsrc
                    Sets the source reference on Chart Studio Cloud
                    for `ticklabelindex`.
                ticklabelmode
                    Determines where tick labels are drawn with
                    respect to their corresponding ticks and grid
                    lines. Only has an effect for axes of `type`
                    "date" When set to "period", tick labels are
                    drawn in the middle of the period between
                    ticks.
                ticklabeloverflow
                    Determines how we handle tick labels that would
                    overflow either the graph div or the domain of
                    the axis. The default value for inside tick
                    labels is *hide past domain*. Otherwise on
                    "category" and "multicategory" axes the default
                    is "allow". In other cases the default is *hide
                    past div*.
                ticklabelposition
                    Determines where tick labels are drawn with
                    respect to the axis Please note that top or
                    bottom has no effect on x axes or when
                    `ticklabelmode` is set to "period". Similarly
                    left or right has no effect on y axes or when
                    `ticklabelmode` is set to "period". Has no
                    effect on "multicategory" axes or when
                    `tickson` is set to "boundaries". When used on
                    axes linked by `matches` or `scaleanchor`, no
                    extra padding for inside labels would be added
                    by autorange, so that the scales could match.
                ticklabelshift
                    Shifts the tick labels by the specified number
                    of pixels in parallel to the axis. Positive
                    values move the labels in the positive
                    direction of the axis.
                ticklabelstandoff
                    Sets the standoff distance (in px) between the
                    axis tick labels and their default position. A
                    positive `ticklabelstandoff` moves the labels
                    farther away from the plot area if
                    `ticklabelposition` is "outside", and deeper
                    into the plot area if `ticklabelposition` is
                    "inside". A negative `ticklabelstandoff` works
                    in the opposite direction, moving outside ticks
                    towards the plot area and inside ticks towards
                    the outside. If the negative value is large
                    enough, inside ticks can even end up outside
                    and vice versa.
                ticklabelstep
                    Sets the spacing between tick labels as
                    compared to the spacing between ticks. A value
                    of 1 (default) means each tick gets a label. A
                    value of 2 means shows every 2nd label. A
                    larger value n means only every nth tick is
                    labeled. `tick0` determines which labels are
                    shown. Not implemented for axes with `type`
                    "log" or "multicategory", or when `tickmode` is
                    "array".
                ticklen
                    Sets the tick length (in px).
                tickmode
                    Sets the tick mode for this axis. If "auto",
                    the number of ticks is set via `nticks`. If
                    "linear", the placement of the ticks is
                    determined by a starting position `tick0` and a
                    tick step `dtick` ("linear" is the default
                    value if `tick0` and `dtick` are provided). If
                    "array", the placement of the ticks is set via
                    `tickvals` and the tick text is `ticktext`.
                    ("array" is the default value if `tickvals` is
                    provided). If "sync", the number of ticks will
                    sync with the overlayed axis set by
                    `overlaying` property.
                tickprefix
                    Sets a tick label prefix.
                ticks
                    Determines whether ticks are drawn or not. If
                    "", this axis' ticks are not drawn. If
                    "outside" ("inside"), this axis' are drawn
                    outside (inside) the axis lines.
                tickson
                    Determines where ticks and grid lines are drawn
                    with respect to their corresponding tick
                    labels. Only has an effect for axes of `type`
                    "category" or "multicategory". When set to
                    "boundaries", ticks and grid lines are drawn
                    half a category to the left/bottom of labels.
                ticksuffix
                    Sets a tick label suffix.
                ticktext
                    Sets the text displayed at the ticks position
                    via `tickvals`. Only has an effect if
                    `tickmode` is set to "array". Used with
                    `tickvals`.
                ticktextsrc
                    Sets the source reference on Chart Studio Cloud
                    for `ticktext`.
                tickvals
                    Sets the values at which ticks on this axis
                    appear. Only has an effect if `tickmode` is set
                    to "array". Used with `ticktext`.
                tickvalssrc
                    Sets the source reference on Chart Studio Cloud
                    for `tickvals`.
                tickwidth
                    Sets the tick width (in px).
                title
                    :class:`plotly.graph_objects.layout.xaxis.Title
                    ` instance or dict with compatible properties
                titlefont
                    Deprecated: Please use layout.xaxis.title.font
                    instead. Sets this axis' title font. Note that
                    the title's font used to be customized by the
                    now deprecated `titlefont` attribute.
                type
                    Sets the axis type. By default, plotly attempts
                    to determined the axis type by looking into the
                    data of the traces that referenced the axis in
                    question.
                uirevision
                    Controls persistence of user-driven changes in
                    axis `range`, `autorange`, and `title` if in
                    `editable: true` configuration. Defaults to
                    `layout.uirevision`.
                visible
                    A single toggle to hide the axis while
                    preserving interaction like dragging. Default
                    is true when a cheater plot is present on the
                    axis, otherwise false
                zeroline
                    Determines whether or not a line is drawn at
                    along the 0 value of this axis. If True, the
                    zero line is drawn on top of the grid lines.
                zerolinecolor
                    Sets the line color of the zero line.
                zerolinewidth
                    Sets the width (in px) of the zero line.

        Returns
        -------
        plotly.graph_objs.layout.XAxis
        """
        return self["xaxis"]

    @xaxis.setter
    def xaxis(self, val):
        self["xaxis"] = val

    # yaxis
    # -----
    @property
    def yaxis(self):
        """
        The 'yaxis' property is an instance of YAxis
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.YAxis`
          - A dict of string/value properties that will be passed
            to the YAxis constructor

            Supported dict properties:

                anchor
                    If set to an opposite-letter axis id (e.g.
                    `x2`, `y`), this axis is bound to the
                    corresponding opposite-letter axis. If set to
                    "free", this axis' position is determined by
                    `position`.
                automargin
                    Determines whether long tick labels
                    automatically grow the figure margins.
                autorange
                    Determines whether or not the range of this
                    axis is computed in relation to the input data.
                    See `rangemode` for more info. If `range` is
                    provided and it has a value for both the lower
                    and upper bound, `autorange` is set to False.
                    Using "min" applies autorange only to set the
                    minimum. Using "max" applies autorange only to
                    set the maximum. Using *min reversed* applies
                    autorange only to set the minimum on a reversed
                    axis. Using *max reversed* applies autorange
                    only to set the maximum on a reversed axis.
                    Using "reversed" applies autorange on both ends
                    and reverses the axis direction.
                autorangeoptions
                    :class:`plotly.graph_objects.layout.yaxis.Autor
                    angeoptions` instance or dict with compatible
                    properties
                autoshift
                    Automatically reposition the axis to avoid
                    overlap with other axes with the same
                    `overlaying` value. This repositioning will
                    account for any `shift` amount applied to other
                    axes on the same side with `autoshift` is set
                    to true. Only has an effect if `anchor` is set
                    to "free".
                autotickangles
                    When `tickangle` is set to "auto", it will be
                    set to the first angle in this array that is
                    large enough to prevent label overlap.
                autotypenumbers
                    Using "strict" a numeric string in trace data
                    is not converted to a number. Using *convert
                    types* a numeric string in trace data may be
                    treated as a number during automatic axis
                    `type` detection. Defaults to
                    layout.autotypenumbers.
                calendar
                    Sets the calendar system to use for `range` and
                    `tick0` if this is a date axis. This does not
                    set the calendar for interpreting data on this
                    axis, that's specified in the trace or via the
                    global `layout.calendar`
                categoryarray
                    Sets the order in which categories on this axis
                    appear. Only has an effect if `categoryorder`
                    is set to "array". Used with `categoryorder`.
                categoryarraysrc
                    Sets the source reference on Chart Studio Cloud
                    for `categoryarray`.
                categoryorder
                    Specifies the ordering logic for the case of
                    categorical variables. By default, plotly uses
                    "trace", which specifies the order that is
                    present in the data supplied. Set
                    `categoryorder` to *category ascending* or
                    *category descending* if order should be
                    determined by the alphanumerical order of the
                    category names. Set `categoryorder` to "array"
                    to derive the ordering from the attribute
                    `categoryarray`. If a category is not found in
                    the `categoryarray` array, the sorting behavior
                    for that attribute will be identical to the
                    "trace" mode. The unspecified categories will
                    follow the categories in `categoryarray`. Set
                    `categoryorder` to *total ascending* or *total
                    descending* if order should be determined by
                    the numerical order of the values. Similarly,
                    the order can be determined by the min, max,
                    sum, mean, geometric mean or median of all the
                    values.
                color
                    Sets default for all colors associated with
                    this axis all at once: line, font, tick, and
                    grid colors. Grid color is lightened by
                    blending this with the plot background
                    Individual pieces can override this.
                constrain
                    If this axis needs to be compressed (either due
                    to its own `scaleanchor` and `scaleratio` or
                    those of the other axis), determines how that
                    happens: by increasing the "range", or by
                    decreasing the "domain". Default is "domain"
                    for axes containing image traces, "range"
                    otherwise.
                constraintoward
                    If this axis needs to be compressed (either due
                    to its own `scaleanchor` and `scaleratio` or
                    those of the other axis), determines which
                    direction we push the originally specified plot
                    area. Options are "left", "center" (default),
                    and "right" for x axes, and "top", "middle"
                    (default), and "bottom" for y axes.
                dividercolor
                    Sets the color of the dividers Only has an
                    effect on "multicategory" axes.
                dividerwidth
                    Sets the width (in px) of the dividers Only has
                    an effect on "multicategory" axes.
                domain
                    Sets the domain of this axis (in plot
                    fraction).
                dtick
                    Sets the step in-between ticks on this axis.
                    Use with `tick0`. Must be a positive number, or
                    special strings available to "log" and "date"
                    axes. If the axis `type` is "log", then ticks
                    are set every 10^(n*dtick) where n is the tick
                    number. For example, to set a tick mark at 1,
                    10, 100, 1000, ... set dtick to 1. To set tick
                    marks at 1, 100, 10000, ... set dtick to 2. To
                    set tick marks at 1, 5, 25, 125, 625, 3125, ...
                    set dtick to log_10(5), or 0.69897000433. "log"
                    has several special values; "L<f>", where `f`
                    is a positive number, gives ticks linearly
                    spaced in value (but not position). For example
                    `tick0` = 0.1, `dtick` = "L0.5" will put ticks
                    at 0.1, 0.6, 1.1, 1.6 etc. To show powers of 10
                    plus small digits between, use "D1" (all
                    digits) or "D2" (only 2 and 5). `tick0` is
                    ignored for "D1" and "D2". If the axis `type`
                    is "date", then you must convert the time to
                    milliseconds. For example, to set the interval
                    between ticks to one day, set `dtick` to
                    86400000.0. "date" also has special values
                    "M<n>" gives ticks spaced by a number of
                    months. `n` must be a positive integer. To set
                    ticks on the 15th of every third month, set
                    `tick0` to "2000-01-15" and `dtick` to "M3". To
                    set ticks every 4 years, set `dtick` to "M48"
                exponentformat
                    Determines a formatting rule for the tick
                    exponents. For example, consider the number
                    1,000,000,000. If "none", it appears as
                    1,000,000,000. If "e", 1e+9. If "E", 1E+9. If
                    "power", 1x10^9 (with 9 in a super script). If
                    "SI", 1G. If "B", 1B.
                fixedrange
                    Determines whether or not this axis is zoom-
                    able. If true, then zoom is disabled.
                gridcolor
                    Sets the color of the grid lines.
                griddash
                    Sets the dash style of lines. Set to a dash
                    type string ("solid", "dot", "dash",
                    "longdash", "dashdot", or "longdashdot") or a
                    dash length list in px (eg "5px,10px,2px,2px").
                gridwidth
                    Sets the width (in px) of the grid lines.
                hoverformat
                    Sets the hover text formatting rule using d3
                    formatting mini-languages which are very
                    similar to those in Python. For numbers, see: h
                    ttps://github.com/d3/d3-format/tree/v1.4.5#d3-
                    format. And for dates see:
                    https://github.com/d3/d3-time-
                    format/tree/v2.2.3#locale_format. We add two
                    items to d3's date formatter: "%h" for half of
                    the year as a decimal number as well as "%{n}f"
                    for fractional seconds with n digits. For
                    example, *2016-10-13 09:15:23.456* with
                    tickformat "%H~%M~%S.%2f" would display
                    "09~15~23.46"
                insiderange
                    Could be used to set the desired inside range
                    of this axis (excluding the labels) when
                    `ticklabelposition` of the anchored axis has
                    "inside". Not implemented for axes with `type`
                    "log". This would be ignored when `range` is
                    provided.
                labelalias
                    Replacement text for specific tick or hover
                    labels. For example using {US: 'USA', CA:
                    'Canada'} changes US to USA and CA to Canada.
                    The labels we would have shown must match the
                    keys exactly, after adding any tickprefix or
                    ticksuffix. For negative numbers the minus sign
                    symbol used (U+2212) is wider than the regular
                    ascii dash. That means you need to use −1
                    instead of -1. labelalias can be used with any
                    axis type, and both keys (if needed) and values
                    (if desired) can include html-like tags or
                    MathJax.
                layer
                    Sets the layer on which this axis is displayed.
                    If *above traces*, this axis is displayed above
                    all the subplot's traces If *below traces*,
                    this axis is displayed below all the subplot's
                    traces, but above the grid lines. Useful when
                    used together with scatter-like traces with
                    `cliponaxis` set to False to show markers
                    and/or text nodes above this axis.
                linecolor
                    Sets the axis line color.
                linewidth
                    Sets the width (in px) of the axis line.
                matches
                    If set to another axis id (e.g. `x2`, `y`), the
                    range of this axis will match the range of the
                    corresponding axis in data-coordinates space.
                    Moreover, matching axes share auto-range
                    values, category lists and histogram auto-bins.
                    Note that setting axes simultaneously in both a
                    `scaleanchor` and a `matches` constraint is
                    currently forbidden. Moreover, note that
                    matching axes must have the same `type`.
                maxallowed
                    Determines the maximum range of this axis.
                minallowed
                    Determines the minimum range of this axis.
                minexponent
                    Hide SI prefix for 10^n if |n| is below this
                    number. This only has an effect when
                    `tickformat` is "SI" or "B".
                minor
                    :class:`plotly.graph_objects.layout.yaxis.Minor
                    ` instance or dict with compatible properties
                mirror
                    Determines if the axis lines or/and ticks are
                    mirrored to the opposite side of the plotting
                    area. If True, the axis lines are mirrored. If
                    "ticks", the axis lines and ticks are mirrored.
                    If False, mirroring is disable. If "all", axis
                    lines are mirrored on all shared-axes subplots.
                    If "allticks", axis lines and ticks are
                    mirrored on all shared-axes subplots.
                nticks
                    Specifies the maximum number of ticks for the
                    particular axis. The actual number of ticks
                    will be chosen automatically to be less than or
                    equal to `nticks`. Has an effect only if
                    `tickmode` is set to "auto".
                overlaying
                    If set a same-letter axis id, this axis is
                    overlaid on top of the corresponding same-
                    letter axis, with traces and axes visible for
                    both axes. If False, this axis does not overlay
                    any same-letter axes. In this case, for axes
                    with overlapping domains only the highest-
                    numbered axis will be visible.
                position
                    Sets the position of this axis in the plotting
                    space (in normalized coordinates). Only has an
                    effect if `anchor` is set to "free".
                range
                    Sets the range of this axis. If the axis `type`
                    is "log", then you must take the log of your
                    desired range (e.g. to set the range from 1 to
                    100, set the range from 0 to 2). If the axis
                    `type` is "date", it should be date strings,
                    like date data, though Date objects and unix
                    milliseconds will be accepted and converted to
                    strings. If the axis `type` is "category", it
                    should be numbers, using the scale where each
                    category is assigned a serial number from zero
                    in the order it appears. Leaving either or both
                    elements `null` impacts the default
                    `autorange`.
                rangebreaks
                    A tuple of :class:`plotly.graph_objects.layout.
                    yaxis.Rangebreak` instances or dicts with
                    compatible properties
                rangebreakdefaults
                    When used in a template (as layout.template.lay
                    out.yaxis.rangebreakdefaults), sets the default
                    property values to use for elements of
                    layout.yaxis.rangebreaks
                rangemode
                    If "normal", the range is computed in relation
                    to the extrema of the input data. If *tozero*`,
                    the range extends to 0, regardless of the input
                    data If "nonnegative", the range is non-
                    negative, regardless of the input data. Applies
                    only to linear axes.
                scaleanchor
                    If set to another axis id (e.g. `x2`, `y`), the
                    range of this axis changes together with the
                    range of the corresponding axis such that the
                    scale of pixels per unit is in a constant
                    ratio. Both axes are still zoomable, but when
                    you zoom one, the other will zoom the same
                    amount, keeping a fixed midpoint. `constrain`
                    and `constraintoward` determine how we enforce
                    the constraint. You can chain these, ie `yaxis:
                    {scaleanchor: *x*}, xaxis2: {scaleanchor: *y*}`
                    but you can only link axes of the same `type`.
                    The linked axis can have the opposite letter
                    (to constrain the aspect ratio) or the same
                    letter (to match scales across subplots). Loops
                    (`yaxis: {scaleanchor: *x*}, xaxis:
                    {scaleanchor: *y*}` or longer) are redundant
                    and the last constraint encountered will be
                    ignored to avoid possible inconsistent
                    constraints via `scaleratio`. Note that setting
                    axes simultaneously in both a `scaleanchor` and
                    a `matches` constraint is currently forbidden.
                    Setting `false` allows to remove a default
                    constraint (occasionally, you may need to
                    prevent a default `scaleanchor` constraint from
                    being applied, eg. when having an image trace
                    `yaxis: {scaleanchor: "x"}` is set
                    automatically in order for pixels to be
                    rendered as squares, setting `yaxis:
                    {scaleanchor: false}` allows to remove the
                    constraint).
                scaleratio
                    If this axis is linked to another by
                    `scaleanchor`, this determines the pixel to
                    unit scale ratio. For example, if this value is
                    10, then every unit on this axis spans 10 times
                    the number of pixels as a unit on the linked
                    axis. Use this for example to create an
                    elevation profile where the vertical scale is
                    exaggerated a fixed amount with respect to the
                    horizontal.
                separatethousands
                    If "true", even 4-digit integers are separated
                shift
                    Moves the axis a given number of pixels from
                    where it would have been otherwise. Accepts
                    both positive and negative values, which will
                    shift the axis either right or left,
                    respectively. If `autoshift` is set to true,
                    then this defaults to a padding of -3 if `side`
                    is set to "left". and defaults to +3 if `side`
                    is set to "right". Defaults to 0 if `autoshift`
                    is set to false. Only has an effect if `anchor`
                    is set to "free".
                showdividers
                    Determines whether or not a dividers are drawn
                    between the category levels of this axis. Only
                    has an effect on "multicategory" axes.
                showexponent
                    If "all", all exponents are shown besides their
                    significands. If "first", only the exponent of
                    the first tick is shown. If "last", only the
                    exponent of the last tick is shown. If "none",
                    no exponents appear.
                showgrid
                    Determines whether or not grid lines are drawn.
                    If True, the grid lines are drawn at every tick
                    mark.
                showline
                    Determines whether or not a line bounding this
                    axis is drawn.
                showspikes
                    Determines whether or not spikes (aka
                    droplines) are drawn for this axis. Note: This
                    only takes affect when hovermode = closest
                showticklabels
                    Determines whether or not the tick labels are
                    drawn.
                showtickprefix
                    If "all", all tick labels are displayed with a
                    prefix. If "first", only the first tick is
                    displayed with a prefix. If "last", only the
                    last tick is displayed with a suffix. If
                    "none", tick prefixes are hidden.
                showticksuffix
                    Same as `showtickprefix` but for tick suffixes.
                side
                    Determines whether a x (y) axis is positioned
                    at the "bottom" ("left") or "top" ("right") of
                    the plotting area.
                spikecolor
                    Sets the spike color. If undefined, will use
                    the series color
                spikedash
                    Sets the dash style of lines. Set to a dash
                    type string ("solid", "dot", "dash",
                    "longdash", "dashdot", or "longdashdot") or a
                    dash length list in px (eg "5px,10px,2px,2px").
                spikemode
                    Determines the drawing mode for the spike line
                    If "toaxis", the line is drawn from the data
                    point to the axis the  series is plotted on. If
                    "across", the line is drawn across the entire
                    plot area, and supercedes "toaxis". If
                    "marker", then a marker dot is drawn on the
                    axis the series is plotted on
                spikesnap
                    Determines whether spikelines are stuck to the
                    cursor or to the closest datapoints.
                spikethickness
                    Sets the width (in px) of the zero line.
                tick0
                    Sets the placement of the first tick on this
                    axis. Use with `dtick`. If the axis `type` is
                    "log", then you must take the log of your
                    starting tick (e.g. to set the starting tick to
                    100, set the `tick0` to 2) except when
                    `dtick`=*L<f>* (see `dtick` for more info). If
                    the axis `type` is "date", it should be a date
                    string, like date data. If the axis `type` is
                    "category", it should be a number, using the
                    scale where each category is assigned a serial
                    number from zero in the order it appears.
                tickangle
                    Sets the angle of the tick labels with respect
                    to the horizontal. For example, a `tickangle`
                    of -90 draws the tick labels vertically.
                tickcolor
                    Sets the tick color.
                tickfont
                    Sets the tick font.
                tickformat
                    Sets the tick label formatting rule using d3
                    formatting mini-languages which are very
                    similar to those in Python. For numbers, see: h
                    ttps://github.com/d3/d3-format/tree/v1.4.5#d3-
                    format. And for dates see:
                    https://github.com/d3/d3-time-
                    format/tree/v2.2.3#locale_format. We add two
                    items to d3's date formatter: "%h" for half of
                    the year as a decimal number as well as "%{n}f"
                    for fractional seconds with n digits. For
                    example, *2016-10-13 09:15:23.456* with
                    tickformat "%H~%M~%S.%2f" would display
                    "09~15~23.46"
                tickformatstops
                    A tuple of :class:`plotly.graph_objects.layout.
                    yaxis.Tickformatstop` instances or dicts with
                    compatible properties
                tickformatstopdefaults
                    When used in a template (as layout.template.lay
                    out.yaxis.tickformatstopdefaults), sets the
                    default property values to use for elements of
                    layout.yaxis.tickformatstops
                ticklabelindex
                    Only for axes with `type` "date" or "linear".
                    Instead of drawing the major tick label, draw
                    the label for the minor tick that is n
                    positions away from the major tick. E.g. to
                    always draw the label for the minor tick before
                    each major tick, choose `ticklabelindex` -1.
                    This is useful for date axes with
                    `ticklabelmode` "period" if you want to label
                    the period that ends with each major tick
                    instead of the period that begins there.
                ticklabelindexsrc
                    Sets the source reference on Chart Studio Cloud
                    for `ticklabelindex`.
                ticklabelmode
                    Determines where tick labels are drawn with
                    respect to their corresponding ticks and grid
                    lines. Only has an effect for axes of `type`
                    "date" When set to "period", tick labels are
                    drawn in the middle of the period between
                    ticks.
                ticklabeloverflow
                    Determines how we handle tick labels that would
                    overflow either the graph div or the domain of
                    the axis. The default value for inside tick
                    labels is *hide past domain*. Otherwise on
                    "category" and "multicategory" axes the default
                    is "allow". In other cases the default is *hide
                    past div*.
                ticklabelposition
                    Determines where tick labels are drawn with
                    respect to the axis Please note that top or
                    bottom has no effect on x axes or when
                    `ticklabelmode` is set to "period". Similarly
                    left or right has no effect on y axes or when
                    `ticklabelmode` is set to "period". Has no
                    effect on "multicategory" axes or when
                    `tickson` is set to "boundaries". When used on
                    axes linked by `matches` or `scaleanchor`, no
                    extra padding for inside labels would be added
                    by autorange, so that the scales could match.
                ticklabelshift
                    Shifts the tick labels by the specified number
                    of pixels in parallel to the axis. Positive
                    values move the labels in the positive
                    direction of the axis.
                ticklabelstandoff
                    Sets the standoff distance (in px) between the
                    axis tick labels and their default position. A
                    positive `ticklabelstandoff` moves the labels
                    farther away from the plot area if
                    `ticklabelposition` is "outside", and deeper
                    into the plot area if `ticklabelposition` is
                    "inside". A negative `ticklabelstandoff` works
                    in the opposite direction, moving outside ticks
                    towards the plot area and inside ticks towards
                    the outside. If the negative value is large
                    enough, inside ticks can even end up outside
                    and vice versa.
                ticklabelstep
                    Sets the spacing between tick labels as
                    compared to the spacing between ticks. A value
                    of 1 (default) means each tick gets a label. A
                    value of 2 means shows every 2nd label. A
                    larger value n means only every nth tick is
                    labeled. `tick0` determines which labels are
                    shown. Not implemented for axes with `type`
                    "log" or "multicategory", or when `tickmode` is
                    "array".
                ticklen
                    Sets the tick length (in px).
                tickmode
                    Sets the tick mode for this axis. If "auto",
                    the number of ticks is set via `nticks`. If
                    "linear", the placement of the ticks is
                    determined by a starting position `tick0` and a
                    tick step `dtick` ("linear" is the default
                    value if `tick0` and `dtick` are provided). If
                    "array", the placement of the ticks is set via
                    `tickvals` and the tick text is `ticktext`.
                    ("array" is the default value if `tickvals` is
                    provided). If "sync", the number of ticks will
                    sync with the overlayed axis set by
                    `overlaying` property.
                tickprefix
                    Sets a tick label prefix.
                ticks
                    Determines whether ticks are drawn or not. If
                    "", this axis' ticks are not drawn. If
                    "outside" ("inside"), this axis' are drawn
                    outside (inside) the axis lines.
                tickson
                    Determines where ticks and grid lines are drawn
                    with respect to their corresponding tick
                    labels. Only has an effect for axes of `type`
                    "category" or "multicategory". When set to
                    "boundaries", ticks and grid lines are drawn
                    half a category to the left/bottom of labels.
                ticksuffix
                    Sets a tick label suffix.
                ticktext
                    Sets the text displayed at the ticks position
                    via `tickvals`. Only has an effect if
                    `tickmode` is set to "array". Used with
                    `tickvals`.
                ticktextsrc
                    Sets the source reference on Chart Studio Cloud
                    for `ticktext`.
                tickvals
                    Sets the values at which ticks on this axis
                    appear. Only has an effect if `tickmode` is set
                    to "array". Used with `ticktext`.
                tickvalssrc
                    Sets the source reference on Chart Studio Cloud
                    for `tickvals`.
                tickwidth
                    Sets the tick width (in px).
                title
                    :class:`plotly.graph_objects.layout.yaxis.Title
                    ` instance or dict with compatible properties
                titlefont
                    Deprecated: Please use layout.yaxis.title.font
                    instead. Sets this axis' title font. Note that
                    the title's font used to be customized by the
                    now deprecated `titlefont` attribute.
                type
                    Sets the axis type. By default, plotly attempts
                    to determined the axis type by looking into the
                    data of the traces that referenced the axis in
                    question.
                uirevision
                    Controls persistence of user-driven changes in
                    axis `range`, `autorange`, and `title` if in
                    `editable: true` configuration. Defaults to
                    `layout.uirevision`.
                visible
                    A single toggle to hide the axis while
                    preserving interaction like dragging. Default
                    is true when a cheater plot is present on the
                    axis, otherwise false
                zeroline
                    Determines whether or not a line is drawn at
                    along the 0 value of this axis. If True, the
                    zero line is drawn on top of the grid lines.
                zerolinecolor
                    Sets the line color of the zero line.
                zerolinewidth
                    Sets the width (in px) of the zero line.

        Returns
        -------
        plotly.graph_objs.layout.YAxis
        """
        return self["yaxis"]

    @yaxis.setter
    def yaxis(self, val):
        self["yaxis"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        activeselection
            :class:`plotly.graph_objects.layout.Activeselection`
            instance or dict with compatible properties
        activeshape
            :class:`plotly.graph_objects.layout.Activeshape`
            instance or dict with compatible properties
        annotations
            A tuple of
            :class:`plotly.graph_objects.layout.Annotation`
            instances or dicts with compatible properties
        annotationdefaults
            When used in a template (as
            layout.template.layout.annotationdefaults), sets the
            default property values to use for elements of
            layout.annotations
        autosize
            Determines whether or not a layout width or height that
            has been left undefined by the user is initialized on
            each relayout. Note that, regardless of this attribute,
            an undefined layout width or height is always
            initialized on the first call to plot.
        autotypenumbers
            Using "strict" a numeric string in trace data is not
            converted to a number. Using *convert types* a numeric
            string in trace data may be treated as a number during
            automatic axis `type` detection. This is the default
            value; however it could be overridden for individual
            axes.
        barcornerradius
            Sets the rounding of bar corners. May be an integer
            number of pixels, or a percentage of bar width (as a
            string ending in %).
        bargap
            Sets the gap (in plot fraction) between bars of
            adjacent location coordinates.
        bargroupgap
            Sets the gap (in plot fraction) between bars of the
            same location coordinate.
        barmode
            Determines how bars at the same location coordinate are
            displayed on the graph. With "stack", the bars are
            stacked on top of one another With "relative", the bars
            are stacked on top of one another, with negative values
            below the axis, positive values above With "group", the
            bars are plotted next to one another centered around
            the shared location. With "overlay", the bars are
            plotted over one another, you might need to reduce
            "opacity" to see multiple bars.
        barnorm
            Sets the normalization for bar traces on the graph.
            With "fraction", the value of each bar is divided by
            the sum of all values at that location coordinate.
            "percent" is the same but multiplied by 100 to show
            percentages.
        boxgap
            Sets the gap (in plot fraction) between boxes of
            adjacent location coordinates. Has no effect on traces
            that have "width" set.
        boxgroupgap
            Sets the gap (in plot fraction) between boxes of the
            same location coordinate. Has no effect on traces that
            have "width" set.
        boxmode
            Determines how boxes at the same location coordinate
            are displayed on the graph. If "group", the boxes are
            plotted next to one another centered around the shared
            location. If "overlay", the boxes are plotted over one
            another, you might need to set "opacity" to see them
            multiple boxes. Has no effect on traces that have
            "width" set.
        calendar
            Sets the default calendar system to use for
            interpreting and displaying dates throughout the plot.
        clickmode
            Determines the mode of single click interactions.
            "event" is the default value and emits the
            `plotly_click` event. In addition this mode emits the
            `plotly_selected` event in drag modes "lasso" and
            "select", but with no event data attached (kept for
            compatibility reasons). The "select" flag enables
            selecting single data points via click. This mode also
            supports persistent selections, meaning that pressing
            Shift while clicking, adds to / subtracts from an
            existing selection. "select" with `hovermode`: "x" can
            be confusing, consider explicitly setting `hovermode`:
            "closest" when using this feature. Selection events are
            sent accordingly as long as "event" flag is set as
            well. When the "event" flag is missing, `plotly_click`
            and `plotly_selected` events are not fired.
        coloraxis
            :class:`plotly.graph_objects.layout.Coloraxis` instance
            or dict with compatible properties
        colorscale
            :class:`plotly.graph_objects.layout.Colorscale`
            instance or dict with compatible properties
        colorway
            Sets the default trace colors.
        computed
            Placeholder for exporting automargin-impacting values
            namely `margin.t`, `margin.b`, `margin.l` and
            `margin.r` in "full-json" mode.
        datarevision
            If provided, a changed value tells `Plotly.react` that
            one or more data arrays has changed. This way you can
            modify arrays in-place rather than making a complete
            new copy for an incremental change. If NOT provided,
            `Plotly.react` assumes that data arrays are being
            treated as immutable, thus any data array with a
            different identity from its predecessor contains new
            data.
        dragmode
            Determines the mode of drag interactions. "select" and
            "lasso" apply only to scatter traces with markers or
            text. "orbit" and "turntable" apply only to 3D scenes.
        editrevision
            Controls persistence of user-driven changes in
            `editable: true` configuration, other than trace names
            and axis titles. Defaults to `layout.uirevision`.
        extendfunnelareacolors
            If `true`, the funnelarea slice colors (whether given
            by `funnelareacolorway` or inherited from `colorway`)
            will be extended to three times its original length by
            first repeating every color 20% lighter then each color
            20% darker. This is intended to reduce the likelihood
            of reusing the same color when you have many slices,
            but you can set `false` to disable. Colors provided in
            the trace, using `marker.colors`, are never extended.
        extendiciclecolors
            If `true`, the icicle slice colors (whether given by
            `iciclecolorway` or inherited from `colorway`) will be
            extended to three times its original length by first
            repeating every color 20% lighter then each color 20%
            darker. This is intended to reduce the likelihood of
            reusing the same color when you have many slices, but
            you can set `false` to disable. Colors provided in the
            trace, using `marker.colors`, are never extended.
        extendpiecolors
            If `true`, the pie slice colors (whether given by
            `piecolorway` or inherited from `colorway`) will be
            extended to three times its original length by first
            repeating every color 20% lighter then each color 20%
            darker. This is intended to reduce the likelihood of
            reusing the same color when you have many slices, but
            you can set `false` to disable. Colors provided in the
            trace, using `marker.colors`, are never extended.
        extendsunburstcolors
            If `true`, the sunburst slice colors (whether given by
            `sunburstcolorway` or inherited from `colorway`) will
            be extended to three times its original length by first
            repeating every color 20% lighter then each color 20%
            darker. This is intended to reduce the likelihood of
            reusing the same color when you have many slices, but
            you can set `false` to disable. Colors provided in the
            trace, using `marker.colors`, are never extended.
        extendtreemapcolors
            If `true`, the treemap slice colors (whether given by
            `treemapcolorway` or inherited from `colorway`) will be
            extended to three times its original length by first
            repeating every color 20% lighter then each color 20%
            darker. This is intended to reduce the likelihood of
            reusing the same color when you have many slices, but
            you can set `false` to disable. Colors provided in the
            trace, using `marker.colors`, are never extended.
        font
            Sets the global font. Note that fonts used in traces
            and other layout components inherit from the global
            font.
        funnelareacolorway
            Sets the default funnelarea slice colors. Defaults to
            the main `colorway` used for trace colors. If you
            specify a new list here it can still be extended with
            lighter and darker colors, see
            `extendfunnelareacolors`.
        funnelgap
            Sets the gap (in plot fraction) between bars of
            adjacent location coordinates.
        funnelgroupgap
            Sets the gap (in plot fraction) between bars of the
            same location coordinate.
        funnelmode
            Determines how bars at the same location coordinate are
            displayed on the graph. With "stack", the bars are
            stacked on top of one another With "group", the bars
            are plotted next to one another centered around the
            shared location. With "overlay", the bars are plotted
            over one another, you might need to reduce "opacity" to
            see multiple bars.
        geo
            :class:`plotly.graph_objects.layout.Geo` instance or
            dict with compatible properties
        grid
            :class:`plotly.graph_objects.layout.Grid` instance or
            dict with compatible properties
        height
            Sets the plot's height (in px).
        hiddenlabels
            hiddenlabels is the funnelarea & pie chart analog of
            visible:'legendonly' but it can contain many labels,
            and can simultaneously hide slices from several
            pies/funnelarea charts
        hiddenlabelssrc
            Sets the source reference on Chart Studio Cloud for
            `hiddenlabels`.
        hidesources
            Determines whether or not a text link citing the data
            source is placed at the bottom-right cored of the
            figure. Has only an effect only on graphs that have
            been generated via forked graphs from the Chart Studio
            Cloud (at https://chart-studio.plotly.com or on-
            premise).
        hoverdistance
            Sets the default distance (in pixels) to look for data
            to add hover labels (-1 means no cutoff, 0 means no
            looking for data). This is only a real distance for
            hovering on point-like objects, like scatter points.
            For area-like objects (bars, scatter fills, etc)
            hovering is on inside the area and off outside, but
            these objects will not supersede hover on point-like
            objects in case of conflict.
        hoverlabel
            :class:`plotly.graph_objects.layout.Hoverlabel`
            instance or dict with compatible properties
        hovermode
            Determines the mode of hover interactions. If
            "closest", a single hoverlabel will appear for the
            "closest" point within the `hoverdistance`. If "x" (or
            "y"), multiple hoverlabels will appear for multiple
            points at the "closest" x- (or y-) coordinate within
            the `hoverdistance`, with the caveat that no more than
            one hoverlabel will appear per trace. If *x unified*
            (or *y unified*), a single hoverlabel will appear
            multiple points at the closest x- (or y-) coordinate
            within the `hoverdistance` with the caveat that no more
            than one hoverlabel will appear per trace. In this
            mode, spikelines are enabled by default perpendicular
            to the specified axis. If false, hover interactions are
            disabled.
        hoversubplots
            Determines expansion of hover effects to other subplots
            If "single" just the axis pair of the primary point is
            included without overlaying subplots. If "overlaying"
            all subplots using the main axis and occupying the same
            space are included. If "axis", also include stacked
            subplots using the same axis when `hovermode` is set to
            "x", *x unified*, "y" or *y unified*.
        iciclecolorway
            Sets the default icicle slice colors. Defaults to the
            main `colorway` used for trace colors. If you specify a
            new list here it can still be extended with lighter and
            darker colors, see `extendiciclecolors`.
        images
            A tuple of :class:`plotly.graph_objects.layout.Image`
            instances or dicts with compatible properties
        imagedefaults
            When used in a template (as
            layout.template.layout.imagedefaults), sets the default
            property values to use for elements of layout.images
        legend
            :class:`plotly.graph_objects.layout.Legend` instance or
            dict with compatible properties
        map
            :class:`plotly.graph_objects.layout.Map` instance or
            dict with compatible properties
        mapbox
            :class:`plotly.graph_objects.layout.Mapbox` instance or
            dict with compatible properties
        margin
            :class:`plotly.graph_objects.layout.Margin` instance or
            dict with compatible properties
        meta
            Assigns extra meta information that can be used in
            various `text` attributes. Attributes such as the
            graph, axis and colorbar `title.text`, annotation
            `text` `trace.name` in legend items, `rangeselector`,
            `updatemenus` and `sliders` `label` text all support
            `meta`. One can access `meta` fields using template
            strings: `%{meta[i]}` where `i` is the index of the
            `meta` item in question. `meta` can also be an object
            for example `{key: value}` which can be accessed
            %{meta[key]}.
        metasrc
            Sets the source reference on Chart Studio Cloud for
            `meta`.
        minreducedheight
            Minimum height of the plot with margin.automargin
            applied (in px)
        minreducedwidth
            Minimum width of the plot with margin.automargin
            applied (in px)
        modebar
            :class:`plotly.graph_objects.layout.Modebar` instance
            or dict with compatible properties
        newselection
            :class:`plotly.graph_objects.layout.Newselection`
            instance or dict with compatible properties
        newshape
            :class:`plotly.graph_objects.layout.Newshape` instance
            or dict with compatible properties
        paper_bgcolor
            Sets the background color of the paper where the graph
            is drawn.
        piecolorway
            Sets the default pie slice colors. Defaults to the main
            `colorway` used for trace colors. If you specify a new
            list here it can still be extended with lighter and
            darker colors, see `extendpiecolors`.
        plot_bgcolor
            Sets the background color of the plotting area in-
            between x and y axes.
        polar
            :class:`plotly.graph_objects.layout.Polar` instance or
            dict with compatible properties
        scattergap
            Sets the gap (in plot fraction) between scatter points
            of adjacent location coordinates. Defaults to `bargap`.
        scattermode
            Determines how scatter points at the same location
            coordinate are displayed on the graph. With "group",
            the scatter points are plotted next to one another
            centered around the shared location. With "overlay",
            the scatter points are plotted over one another, you
            might need to reduce "opacity" to see multiple scatter
            points.
        scene
            :class:`plotly.graph_objects.layout.Scene` instance or
            dict with compatible properties
        selectdirection
            When `dragmode` is set to "select", this limits the
            selection of the drag to horizontal, vertical or
            diagonal. "h" only allows horizontal selection, "v"
            only vertical, "d" only diagonal and "any" sets no
            limit.
        selectionrevision
            Controls persistence of user-driven changes in selected
            points from all traces.
        selections
            A tuple of
            :class:`plotly.graph_objects.layout.Selection`
            instances or dicts with compatible properties
        selectiondefaults
            When used in a template (as
            layout.template.layout.selectiondefaults), sets the
            default property values to use for elements of
            layout.selections
        separators
            Sets the decimal and thousand separators. For example,
            *. * puts a '.' before decimals and a space between
            thousands. In English locales, dflt is ".," but other
            locales may alter this default.
        shapes
            A tuple of :class:`plotly.graph_objects.layout.Shape`
            instances or dicts with compatible properties
        shapedefaults
            When used in a template (as
            layout.template.layout.shapedefaults), sets the default
            property values to use for elements of layout.shapes
        showlegend
            Determines whether or not a legend is drawn. Default is
            `true` if there is a trace to show and any of these: a)
            Two or more traces would by default be shown in the
            legend. b) One pie trace is shown in the legend. c) One
            trace is explicitly given with `showlegend: true`.
        sliders
            A tuple of :class:`plotly.graph_objects.layout.Slider`
            instances or dicts with compatible properties
        sliderdefaults
            When used in a template (as
            layout.template.layout.sliderdefaults), sets the
            default property values to use for elements of
            layout.sliders
        smith
            :class:`plotly.graph_objects.layout.Smith` instance or
            dict with compatible properties
        spikedistance
            Sets the default distance (in pixels) to look for data
            to draw spikelines to (-1 means no cutoff, 0 means no
            looking for data). As with hoverdistance, distance does
            not apply to area-like objects. In addition, some
            objects can be hovered on but will not generate
            spikelines, such as scatter fills.
        sunburstcolorway
            Sets the default sunburst slice colors. Defaults to the
            main `colorway` used for trace colors. If you specify a
            new list here it can still be extended with lighter and
            darker colors, see `extendsunburstcolors`.
        template
            Default attributes to be applied to the plot. This
            should be a dict with format: `{'layout':
            layoutTemplate, 'data': {trace_type: [traceTemplate,
            ...], ...}}` where `layoutTemplate` is a dict matching
            the structure of `figure.layout` and `traceTemplate` is
            a dict matching the structure of the trace with type
            `trace_type` (e.g. 'scatter'). Alternatively, this may
            be specified as an instance of
            plotly.graph_objs.layout.Template.  Trace templates are
            applied cyclically to traces of each type. Container
            arrays (eg `annotations`) have special handling: An
            object ending in `defaults` (eg `annotationdefaults`)
            is applied to each array item. But if an item has a
            `templateitemname` key we look in the template array
            for an item with matching `name` and apply that
            instead. If no matching `name` is found we mark the
            item invisible. Any named template item not referenced
            is appended to the end of the array, so this can be
            used to add a watermark annotation or a logo image, for
            example. To omit one of these items on the plot, make
            an item with matching `templateitemname` and `visible:
            false`.
        ternary
            :class:`plotly.graph_objects.layout.Ternary` instance
            or dict with compatible properties
        title
            :class:`plotly.graph_objects.layout.Title` instance or
            dict with compatible properties
        titlefont
            Deprecated: Please use layout.title.font instead. Sets
            the title font. Note that the title's font used to be
            customized by the now deprecated `titlefont` attribute.
        transition
            Sets transition options used during Plotly.react
            updates.
        treemapcolorway
            Sets the default treemap slice colors. Defaults to the
            main `colorway` used for trace colors. If you specify a
            new list here it can still be extended with lighter and
            darker colors, see `extendtreemapcolors`.
        uirevision
            Used to allow user interactions with the plot to
            persist after `Plotly.react` calls that are unaware of
            these interactions. If `uirevision` is omitted, or if
            it is given and it changed from the previous
            `Plotly.react` call, the exact new figure is used. If
            `uirevision` is truthy and did NOT change, any
            attribute that has been affected by user interactions
            and did not receive a different value in the new figure
            will keep the interaction value. `layout.uirevision`
            attribute serves as the default for `uirevision`
            attributes in various sub-containers. For finer control
            you can set these sub-attributes directly. For example,
            if your app separately controls the data on the x and y
            axes you might set `xaxis.uirevision=*time*` and
            `yaxis.uirevision=*cost*`. Then if only the y data is
            changed, you can update `yaxis.uirevision=*quantity*`
            and the y axis range will reset but the x axis range
            will retain any user-driven zoom.
        uniformtext
            :class:`plotly.graph_objects.layout.Uniformtext`
            instance or dict with compatible properties
        updatemenus
            A tuple of
            :class:`plotly.graph_objects.layout.Updatemenu`
            instances or dicts with compatible properties
        updatemenudefaults
            When used in a template (as
            layout.template.layout.updatemenudefaults), sets the
            default property values to use for elements of
            layout.updatemenus
        violingap
            Sets the gap (in plot fraction) between violins of
            adjacent location coordinates. Has no effect on traces
            that have "width" set.
        violingroupgap
            Sets the gap (in plot fraction) between violins of the
            same location coordinate. Has no effect on traces that
            have "width" set.
        violinmode
            Determines how violins at the same location coordinate
            are displayed on the graph. If "group", the violins are
            plotted next to one another centered around the shared
            location. If "overlay", the violins are plotted over
            one another, you might need to set "opacity" to see
            them multiple violins. Has no effect on traces that
            have "width" set.
        waterfallgap
            Sets the gap (in plot fraction) between bars of
            adjacent location coordinates.
        waterfallgroupgap
            Sets the gap (in plot fraction) between bars of the
            same location coordinate.
        waterfallmode
            Determines how bars at the same location coordinate are
            displayed on the graph. With "group", the bars are
            plotted next to one another centered around the shared
            location. With "overlay", the bars are plotted over one
            another, you might need to reduce "opacity" to see
            multiple bars.
        width
            Sets the plot's width (in px).
        xaxis
            :class:`plotly.graph_objects.layout.XAxis` instance or
            dict with compatible properties
        yaxis
            :class:`plotly.graph_objects.layout.YAxis` instance or
            dict with compatible properties
        """

    _mapped_properties = {"titlefont": ("title", "font")}

    def __init__(
        self,
        arg=None,
        activeselection=None,
        activeshape=None,
        annotations=None,
        annotationdefaults=None,
        autosize=None,
        autotypenumbers=None,
        barcornerradius=None,
        bargap=None,
        bargroupgap=None,
        barmode=None,
        barnorm=None,
        boxgap=None,
        boxgroupgap=None,
        boxmode=None,
        calendar=None,
        clickmode=None,
        coloraxis=None,
        colorscale=None,
        colorway=None,
        computed=None,
        datarevision=None,
        dragmode=None,
        editrevision=None,
        extendfunnelareacolors=None,
        extendiciclecolors=None,
        extendpiecolors=None,
        extendsunburstcolors=None,
        extendtreemapcolors=None,
        font=None,
        funnelareacolorway=None,
        funnelgap=None,
        funnelgroupgap=None,
        funnelmode=None,
        geo=None,
        grid=None,
        height=None,
        hiddenlabels=None,
        hiddenlabelssrc=None,
        hidesources=None,
        hoverdistance=None,
        hoverlabel=None,
        hovermode=None,
        hoversubplots=None,
        iciclecolorway=None,
        images=None,
        imagedefaults=None,
        legend=None,
        map=None,
        mapbox=None,
        margin=None,
        meta=None,
        metasrc=None,
        minreducedheight=None,
        minreducedwidth=None,
        modebar=None,
        newselection=None,
        newshape=None,
        paper_bgcolor=None,
        piecolorway=None,
        plot_bgcolor=None,
        polar=None,
        scattergap=None,
        scattermode=None,
        scene=None,
        selectdirection=None,
        selectionrevision=None,
        selections=None,
        selectiondefaults=None,
        separators=None,
        shapes=None,
        shapedefaults=None,
        showlegend=None,
        sliders=None,
        sliderdefaults=None,
        smith=None,
        spikedistance=None,
        sunburstcolorway=None,
        template=None,
        ternary=None,
        title=None,
        titlefont=None,
        transition=None,
        treemapcolorway=None,
        uirevision=None,
        uniformtext=None,
        updatemenus=None,
        updatemenudefaults=None,
        violingap=None,
        violingroupgap=None,
        violinmode=None,
        waterfallgap=None,
        waterfallgroupgap=None,
        waterfallmode=None,
        width=None,
        xaxis=None,
        yaxis=None,
        **kwargs,
    ):
        """
        Construct a new Layout object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of :class:`plotly.graph_objs.Layout`
        activeselection
            :class:`plotly.graph_objects.layout.Activeselection`
            instance or dict with compatible properties
        activeshape
            :class:`plotly.graph_objects.layout.Activeshape`
            instance or dict with compatible properties
        annotations
            A tuple of
            :class:`plotly.graph_objects.layout.Annotation`
            instances or dicts with compatible properties
        annotationdefaults
            When used in a template (as
            layout.template.layout.annotationdefaults), sets the
            default property values to use for elements of
            layout.annotations
        autosize
            Determines whether or not a layout width or height that
            has been left undefined by the user is initialized on
            each relayout. Note that, regardless of this attribute,
            an undefined layout width or height is always
            initialized on the first call to plot.
        autotypenumbers
            Using "strict" a numeric string in trace data is not
            converted to a number. Using *convert types* a numeric
            string in trace data may be treated as a number during
            automatic axis `type` detection. This is the default
            value; however it could be overridden for individual
            axes.
        barcornerradius
            Sets the rounding of bar corners. May be an integer
            number of pixels, or a percentage of bar width (as a
            string ending in %).
        bargap
            Sets the gap (in plot fraction) between bars of
            adjacent location coordinates.
        bargroupgap
            Sets the gap (in plot fraction) between bars of the
            same location coordinate.
        barmode
            Determines how bars at the same location coordinate are
            displayed on the graph. With "stack", the bars are
            stacked on top of one another With "relative", the bars
            are stacked on top of one another, with negative values
            below the axis, positive values above With "group", the
            bars are plotted next to one another centered around
            the shared location. With "overlay", the bars are
            plotted over one another, you might need to reduce
            "opacity" to see multiple bars.
        barnorm
            Sets the normalization for bar traces on the graph.
            With "fraction", the value of each bar is divided by
            the sum of all values at that location coordinate.
            "percent" is the same but multiplied by 100 to show
            percentages.
        boxgap
            Sets the gap (in plot fraction) between boxes of
            adjacent location coordinates. Has no effect on traces
            that have "width" set.
        boxgroupgap
            Sets the gap (in plot fraction) between boxes of the
            same location coordinate. Has no effect on traces that
            have "width" set.
        boxmode
            Determines how boxes at the same location coordinate
            are displayed on the graph. If "group", the boxes are
            plotted next to one another centered around the shared
            location. If "overlay", the boxes are plotted over one
            another, you might need to set "opacity" to see them
            multiple boxes. Has no effect on traces that have
            "width" set.
        calendar
            Sets the default calendar system to use for
            interpreting and displaying dates throughout the plot.
        clickmode
            Determines the mode of single click interactions.
            "event" is the default value and emits the
            `plotly_click` event. In addition this mode emits the
            `plotly_selected` event in drag modes "lasso" and
            "select", but with no event data attached (kept for
            compatibility reasons). The "select" flag enables
            selecting single data points via click. This mode also
            supports persistent selections, meaning that pressing
            Shift while clicking, adds to / subtracts from an
            existing selection. "select" with `hovermode`: "x" can
            be confusing, consider explicitly setting `hovermode`:
            "closest" when using this feature. Selection events are
            sent accordingly as long as "event" flag is set as
            well. When the "event" flag is missing, `plotly_click`
            and `plotly_selected` events are not fired.
        coloraxis
            :class:`plotly.graph_objects.layout.Coloraxis` instance
            or dict with compatible properties
        colorscale
            :class:`plotly.graph_objects.layout.Colorscale`
            instance or dict with compatible properties
        colorway
            Sets the default trace colors.
        computed
            Placeholder for exporting automargin-impacting values
            namely `margin.t`, `margin.b`, `margin.l` and
            `margin.r` in "full-json" mode.
        datarevision
            If provided, a changed value tells `Plotly.react` that
            one or more data arrays has changed. This way you can
            modify arrays in-place rather than making a complete
            new copy for an incremental change. If NOT provided,
            `Plotly.react` assumes that data arrays are being
            treated as immutable, thus any data array with a
            different identity from its predecessor contains new
            data.
        dragmode
            Determines the mode of drag interactions. "select" and
            "lasso" apply only to scatter traces with markers or
            text. "orbit" and "turntable" apply only to 3D scenes.
        editrevision
            Controls persistence of user-driven changes in
            `editable: true` configuration, other than trace names
            and axis titles. Defaults to `layout.uirevision`.
        extendfunnelareacolors
            If `true`, the funnelarea slice colors (whether given
            by `funnelareacolorway` or inherited from `colorway`)
            will be extended to three times its original length by
            first repeating every color 20% lighter then each color
            20% darker. This is intended to reduce the likelihood
            of reusing the same color when you have many slices,
            but you can set `false` to disable. Colors provided in
            the trace, using `marker.colors`, are never extended.
        extendiciclecolors
            If `true`, the icicle slice colors (whether given by
            `iciclecolorway` or inherited from `colorway`) will be
            extended to three times its original length by first
            repeating every color 20% lighter then each color 20%
            darker. This is intended to reduce the likelihood of
            reusing the same color when you have many slices, but
            you can set `false` to disable. Colors provided in the
            trace, using `marker.colors`, are never extended.
        extendpiecolors
            If `true`, the pie slice colors (whether given by
            `piecolorway` or inherited from `colorway`) will be
            extended to three times its original length by first
            repeating every color 20% lighter then each color 20%
            darker. This is intended to reduce the likelihood of
            reusing the same color when you have many slices, but
            you can set `false` to disable. Colors provided in the
            trace, using `marker.colors`, are never extended.
        extendsunburstcolors
            If `true`, the sunburst slice colors (whether given by
            `sunburstcolorway` or inherited from `colorway`) will
            be extended to three times its original length by first
            repeating every color 20% lighter then each color 20%
            darker. This is intended to reduce the likelihood of
            reusing the same color when you have many slices, but
            you can set `false` to disable. Colors provided in the
            trace, using `marker.colors`, are never extended.
        extendtreemapcolors
            If `true`, the treemap slice colors (whether given by
            `treemapcolorway` or inherited from `colorway`) will be
            extended to three times its original length by first
            repeating every color 20% lighter then each color 20%
            darker. This is intended to reduce the likelihood of
            reusing the same color when you have many slices, but
            you can set `false` to disable. Colors provided in the
            trace, using `marker.colors`, are never extended.
        font
            Sets the global font. Note that fonts used in traces
            and other layout components inherit from the global
            font.
        funnelareacolorway
            Sets the default funnelarea slice colors. Defaults to
            the main `colorway` used for trace colors. If you
            specify a new list here it can still be extended with
            lighter and darker colors, see
            `extendfunnelareacolors`.
        funnelgap
            Sets the gap (in plot fraction) between bars of
            adjacent location coordinates.
        funnelgroupgap
            Sets the gap (in plot fraction) between bars of the
            same location coordinate.
        funnelmode
            Determines how bars at the same location coordinate are
            displayed on the graph. With "stack", the bars are
            stacked on top of one another With "group", the bars
            are plotted next to one another centered around the
            shared location. With "overlay", the bars are plotted
            over one another, you might need to reduce "opacity" to
            see multiple bars.
        geo
            :class:`plotly.graph_objects.layout.Geo` instance or
            dict with compatible properties
        grid
            :class:`plotly.graph_objects.layout.Grid` instance or
            dict with compatible properties
        height
            Sets the plot's height (in px).
        hiddenlabels
            hiddenlabels is the funnelarea & pie chart analog of
            visible:'legendonly' but it can contain many labels,
            and can simultaneously hide slices from several
            pies/funnelarea charts
        hiddenlabelssrc
            Sets the source reference on Chart Studio Cloud for
            `hiddenlabels`.
        hidesources
            Determines whether or not a text link citing the data
            source is placed at the bottom-right cored of the
            figure. Has only an effect only on graphs that have
            been generated via forked graphs from the Chart Studio
            Cloud (at https://chart-studio.plotly.com or on-
            premise).
        hoverdistance
            Sets the default distance (in pixels) to look for data
            to add hover labels (-1 means no cutoff, 0 means no
            looking for data). This is only a real distance for
            hovering on point-like objects, like scatter points.
            For area-like objects (bars, scatter fills, etc)
            hovering is on inside the area and off outside, but
            these objects will not supersede hover on point-like
            objects in case of conflict.
        hoverlabel
            :class:`plotly.graph_objects.layout.Hoverlabel`
            instance or dict with compatible properties
        hovermode
            Determines the mode of hover interactions. If
            "closest", a single hoverlabel will appear for the
            "closest" point within the `hoverdistance`. If "x" (or
            "y"), multiple hoverlabels will appear for multiple
            points at the "closest" x- (or y-) coordinate within
            the `hoverdistance`, with the caveat that no more than
            one hoverlabel will appear per trace. If *x unified*
            (or *y unified*), a single hoverlabel will appear
            multiple points at the closest x- (or y-) coordinate
            within the `hoverdistance` with the caveat that no more
            than one hoverlabel will appear per trace. In this
            mode, spikelines are enabled by default perpendicular
            to the specified axis. If false, hover interactions are
            disabled.
        hoversubplots
            Determines expansion of hover effects to other subplots
            If "single" just the axis pair of the primary point is
            included without overlaying subplots. If "overlaying"
            all subplots using the main axis and occupying the same
            space are included. If "axis", also include stacked
            subplots using the same axis when `hovermode` is set to
            "x", *x unified*, "y" or *y unified*.
        iciclecolorway
            Sets the default icicle slice colors. Defaults to the
            main `colorway` used for trace colors. If you specify a
            new list here it can still be extended with lighter and
            darker colors, see `extendiciclecolors`.
        images
            A tuple of :class:`plotly.graph_objects.layout.Image`
            instances or dicts with compatible properties
        imagedefaults
            When used in a template (as
            layout.template.layout.imagedefaults), sets the default
            property values to use for elements of layout.images
        legend
            :class:`plotly.graph_objects.layout.Legend` instance or
            dict with compatible properties
        map
            :class:`plotly.graph_objects.layout.Map` instance or
            dict with compatible properties
        mapbox
            :class:`plotly.graph_objects.layout.Mapbox` instance or
            dict with compatible properties
        margin
            :class:`plotly.graph_objects.layout.Margin` instance or
            dict with compatible properties
        meta
            Assigns extra meta information that can be used in
            various `text` attributes. Attributes such as the
            graph, axis and colorbar `title.text`, annotation
            `text` `trace.name` in legend items, `rangeselector`,
            `updatemenus` and `sliders` `label` text all support
            `meta`. One can access `meta` fields using template
            strings: `%{meta[i]}` where `i` is the index of the
            `meta` item in question. `meta` can also be an object
            for example `{key: value}` which can be accessed
            %{meta[key]}.
        metasrc
            Sets the source reference on Chart Studio Cloud for
            `meta`.
        minreducedheight
            Minimum height of the plot with margin.automargin
            applied (in px)
        minreducedwidth
            Minimum width of the plot with margin.automargin
            applied (in px)
        modebar
            :class:`plotly.graph_objects.layout.Modebar` instance
            or dict with compatible properties
        newselection
            :class:`plotly.graph_objects.layout.Newselection`
            instance or dict with compatible properties
        newshape
            :class:`plotly.graph_objects.layout.Newshape` instance
            or dict with compatible properties
        paper_bgcolor
            Sets the background color of the paper where the graph
            is drawn.
        piecolorway
            Sets the default pie slice colors. Defaults to the main
            `colorway` used for trace colors. If you specify a new
            list here it can still be extended with lighter and
            darker colors, see `extendpiecolors`.
        plot_bgcolor
            Sets the background color of the plotting area in-
            between x and y axes.
        polar
            :class:`plotly.graph_objects.layout.Polar` instance or
            dict with compatible properties
        scattergap
            Sets the gap (in plot fraction) between scatter points
            of adjacent location coordinates. Defaults to `bargap`.
        scattermode
            Determines how scatter points at the same location
            coordinate are displayed on the graph. With "group",
            the scatter points are plotted next to one another
            centered around the shared location. With "overlay",
            the scatter points are plotted over one another, you
            might need to reduce "opacity" to see multiple scatter
            points.
        scene
            :class:`plotly.graph_objects.layout.Scene` instance or
            dict with compatible properties
        selectdirection
            When `dragmode` is set to "select", this limits the
            selection of the drag to horizontal, vertical or
            diagonal. "h" only allows horizontal selection, "v"
            only vertical, "d" only diagonal and "any" sets no
            limit.
        selectionrevision
            Controls persistence of user-driven changes in selected
            points from all traces.
        selections
            A tuple of
            :class:`plotly.graph_objects.layout.Selection`
            instances or dicts with compatible properties
        selectiondefaults
            When used in a template (as
            layout.template.layout.selectiondefaults), sets the
            default property values to use for elements of
            layout.selections
        separators
            Sets the decimal and thousand separators. For example,
            *. * puts a '.' before decimals and a space between
            thousands. In English locales, dflt is ".," but other
            locales may alter this default.
        shapes
            A tuple of :class:`plotly.graph_objects.layout.Shape`
            instances or dicts with compatible properties
        shapedefaults
            When used in a template (as
            layout.template.layout.shapedefaults), sets the default
            property values to use for elements of layout.shapes
        showlegend
            Determines whether or not a legend is drawn. Default is
            `true` if there is a trace to show and any of these: a)
            Two or more traces would by default be shown in the
            legend. b) One pie trace is shown in the legend. c) One
            trace is explicitly given with `showlegend: true`.
        sliders
            A tuple of :class:`plotly.graph_objects.layout.Slider`
            instances or dicts with compatible properties
        sliderdefaults
            When used in a template (as
            layout.template.layout.sliderdefaults), sets the
            default property values to use for elements of
            layout.sliders
        smith
            :class:`plotly.graph_objects.layout.Smith` instance or
            dict with compatible properties
        spikedistance
            Sets the default distance (in pixels) to look for data
            to draw spikelines to (-1 means no cutoff, 0 means no
            looking for data). As with hoverdistance, distance does
            not apply to area-like objects. In addition, some
            objects can be hovered on but will not generate
            spikelines, such as scatter fills.
        sunburstcolorway
            Sets the default sunburst slice colors. Defaults to the
            main `colorway` used for trace colors. If you specify a
            new list here it can still be extended with lighter and
            darker colors, see `extendsunburstcolors`.
        template
            Default attributes to be applied to the plot. This
            should be a dict with format: `{'layout':
            layoutTemplate, 'data': {trace_type: [traceTemplate,
            ...], ...}}` where `layoutTemplate` is a dict matching
            the structure of `figure.layout` and `traceTemplate` is
            a dict matching the structure of the trace with type
            `trace_type` (e.g. 'scatter'). Alternatively, this may
            be specified as an instance of
            plotly.graph_objs.layout.Template.  Trace templates are
            applied cyclically to traces of each type. Container
            arrays (eg `annotations`) have special handling: An
            object ending in `defaults` (eg `annotationdefaults`)
            is applied to each array item. But if an item has a
            `templateitemname` key we look in the template array
            for an item with matching `name` and apply that
            instead. If no matching `name` is found we mark the
            item invisible. Any named template item not referenced
            is appended to the end of the array, so this can be
            used to add a watermark annotation or a logo image, for
            example. To omit one of these items on the plot, make
            an item with matching `templateitemname` and `visible:
            false`.
        ternary
            :class:`plotly.graph_objects.layout.Ternary` instance
            or dict with compatible properties
        title
            :class:`plotly.graph_objects.layout.Title` instance or
            dict with compatible properties
        titlefont
            Deprecated: Please use layout.title.font instead. Sets
            the title font. Note that the title's font used to be
            customized by the now deprecated `titlefont` attribute.
        transition
            Sets transition options used during Plotly.react
            updates.
        treemapcolorway
            Sets the default treemap slice colors. Defaults to the
            main `colorway` used for trace colors. If you specify a
            new list here it can still be extended with lighter and
            darker colors, see `extendtreemapcolors`.
        uirevision
            Used to allow user interactions with the plot to
            persist after `Plotly.react` calls that are unaware of
            these interactions. If `uirevision` is omitted, or if
            it is given and it changed from the previous
            `Plotly.react` call, the exact new figure is used. If
            `uirevision` is truthy and did NOT change, any
            attribute that has been affected by user interactions
            and did not receive a different value in the new figure
            will keep the interaction value. `layout.uirevision`
            attribute serves as the default for `uirevision`
            attributes in various sub-containers. For finer control
            you can set these sub-attributes directly. For example,
            if your app separately controls the data on the x and y
            axes you might set `xaxis.uirevision=*time*` and
            `yaxis.uirevision=*cost*`. Then if only the y data is
            changed, you can update `yaxis.uirevision=*quantity*`
            and the y axis range will reset but the x axis range
            will retain any user-driven zoom.
        uniformtext
            :class:`plotly.graph_objects.layout.Uniformtext`
            instance or dict with compatible properties
        updatemenus
            A tuple of
            :class:`plotly.graph_objects.layout.Updatemenu`
            instances or dicts with compatible properties
        updatemenudefaults
            When used in a template (as
            layout.template.layout.updatemenudefaults), sets the
            default property values to use for elements of
            layout.updatemenus
        violingap
            Sets the gap (in plot fraction) between violins of
            adjacent location coordinates. Has no effect on traces
            that have "width" set.
        violingroupgap
            Sets the gap (in plot fraction) between violins of the
            same location coordinate. Has no effect on traces that
            have "width" set.
        violinmode
            Determines how violins at the same location coordinate
            are displayed on the graph. If "group", the violins are
            plotted next to one another centered around the shared
            location. If "overlay", the violins are plotted over
            one another, you might need to set "opacity" to see
            them multiple violins. Has no effect on traces that
            have "width" set.
        waterfallgap
            Sets the gap (in plot fraction) between bars of
            adjacent location coordinates.
        waterfallgroupgap
            Sets the gap (in plot fraction) between bars of the
            same location coordinate.
        waterfallmode
            Determines how bars at the same location coordinate are
            displayed on the graph. With "group", the bars are
            plotted next to one another centered around the shared
            location. With "overlay", the bars are plotted over one
            another, you might need to reduce "opacity" to see
            multiple bars.
        width
            Sets the plot's width (in px).
        xaxis
            :class:`plotly.graph_objects.layout.XAxis` instance or
            dict with compatible properties
        yaxis
            :class:`plotly.graph_objects.layout.YAxis` instance or
            dict with compatible properties

        Returns
        -------
        Layout
        """
        super(Layout, self).__init__("layout")

        if "_parent" in kwargs:
            self._parent = kwargs["_parent"]
            return

        # Override _valid_props for instance so that instance can mutate set
        # to support subplot properties (e.g. xaxis2)
        self._valid_props = {
            "activeselection",
            "activeshape",
            "annotationdefaults",
            "annotations",
            "autosize",
            "autotypenumbers",
            "barcornerradius",
            "bargap",
            "bargroupgap",
            "barmode",
            "barnorm",
            "boxgap",
            "boxgroupgap",
            "boxmode",
            "calendar",
            "clickmode",
            "coloraxis",
            "colorscale",
            "colorway",
            "computed",
            "datarevision",
            "dragmode",
            "editrevision",
            "extendfunnelareacolors",
            "extendiciclecolors",
            "extendpiecolors",
            "extendsunburstcolors",
            "extendtreemapcolors",
            "font",
            "funnelareacolorway",
            "funnelgap",
            "funnelgroupgap",
            "funnelmode",
            "geo",
            "grid",
            "height",
            "hiddenlabels",
            "hiddenlabelssrc",
            "hidesources",
            "hoverdistance",
            "hoverlabel",
            "hovermode",
            "hoversubplots",
            "iciclecolorway",
            "imagedefaults",
            "images",
            "legend",
            "map",
            "mapbox",
            "margin",
            "meta",
            "metasrc",
            "minreducedheight",
            "minreducedwidth",
            "modebar",
            "newselection",
            "newshape",
            "paper_bgcolor",
            "piecolorway",
            "plot_bgcolor",
            "polar",
            "scattergap",
            "scattermode",
            "scene",
            "selectdirection",
            "selectiondefaults",
            "selectionrevision",
            "selections",
            "separators",
            "shapedefaults",
            "shapes",
            "showlegend",
            "sliderdefaults",
            "sliders",
            "smith",
            "spikedistance",
            "sunburstcolorway",
            "template",
            "ternary",
            "title",
            "titlefont",
            "transition",
            "treemapcolorway",
            "uirevision",
            "uniformtext",
            "updatemenudefaults",
            "updatemenus",
            "violingap",
            "violingroupgap",
            "violinmode",
            "waterfallgap",
            "waterfallgroupgap",
            "waterfallmode",
            "width",
            "xaxis",
            "yaxis",
        }

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
The first argument to the plotly.graph_objs.Layout
constructor must be a dict or
an instance of :class:`plotly.graph_objs.Layout`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("activeselection", None)
        _v = activeselection if activeselection is not None else _v
        if _v is not None:
            self["activeselection"] = _v
        _v = arg.pop("activeshape", None)
        _v = activeshape if activeshape is not None else _v
        if _v is not None:
            self["activeshape"] = _v
        _v = arg.pop("annotations", None)
        _v = annotations if annotations is not None else _v
        if _v is not None:
            self["annotations"] = _v
        _v = arg.pop("annotationdefaults", None)
        _v = annotationdefaults if annotationdefaults is not None else _v
        if _v is not None:
            self["annotationdefaults"] = _v
        _v = arg.pop("autosize", None)
        _v = autosize if autosize is not None else _v
        if _v is not None:
            self["autosize"] = _v
        _v = arg.pop("autotypenumbers", None)
        _v = autotypenumbers if autotypenumbers is not None else _v
        if _v is not None:
            self["autotypenumbers"] = _v
        _v = arg.pop("barcornerradius", None)
        _v = barcornerradius if barcornerradius is not None else _v
        if _v is not None:
            self["barcornerradius"] = _v
        _v = arg.pop("bargap", None)
        _v = bargap if bargap is not None else _v
        if _v is not None:
            self["bargap"] = _v
        _v = arg.pop("bargroupgap", None)
        _v = bargroupgap if bargroupgap is not None else _v
        if _v is not None:
            self["bargroupgap"] = _v
        _v = arg.pop("barmode", None)
        _v = barmode if barmode is not None else _v
        if _v is not None:
            self["barmode"] = _v
        _v = arg.pop("barnorm", None)
        _v = barnorm if barnorm is not None else _v
        if _v is not None:
            self["barnorm"] = _v
        _v = arg.pop("boxgap", None)
        _v = boxgap if boxgap is not None else _v
        if _v is not None:
            self["boxgap"] = _v
        _v = arg.pop("boxgroupgap", None)
        _v = boxgroupgap if boxgroupgap is not None else _v
        if _v is not None:
            self["boxgroupgap"] = _v
        _v = arg.pop("boxmode", None)
        _v = boxmode if boxmode is not None else _v
        if _v is not None:
            self["boxmode"] = _v
        _v = arg.pop("calendar", None)
        _v = calendar if calendar is not None else _v
        if _v is not None:
            self["calendar"] = _v
        _v = arg.pop("clickmode", None)
        _v = clickmode if clickmode is not None else _v
        if _v is not None:
            self["clickmode"] = _v
        _v = arg.pop("coloraxis", None)
        _v = coloraxis if coloraxis is not None else _v
        if _v is not None:
            self["coloraxis"] = _v
        _v = arg.pop("colorscale", None)
        _v = colorscale if colorscale is not None else _v
        if _v is not None:
            self["colorscale"] = _v
        _v = arg.pop("colorway", None)
        _v = colorway if colorway is not None else _v
        if _v is not None:
            self["colorway"] = _v
        _v = arg.pop("computed", None)
        _v = computed if computed is not None else _v
        if _v is not None:
            self["computed"] = _v
        _v = arg.pop("datarevision", None)
        _v = datarevision if datarevision is not None else _v
        if _v is not None:
            self["datarevision"] = _v
        _v = arg.pop("dragmode", None)
        _v = dragmode if dragmode is not None else _v
        if _v is not None:
            self["dragmode"] = _v
        _v = arg.pop("editrevision", None)
        _v = editrevision if editrevision is not None else _v
        if _v is not None:
            self["editrevision"] = _v
        _v = arg.pop("extendfunnelareacolors", None)
        _v = extendfunnelareacolors if extendfunnelareacolors is not None else _v
        if _v is not None:
            self["extendfunnelareacolors"] = _v
        _v = arg.pop("extendiciclecolors", None)
        _v = extendiciclecolors if extendiciclecolors is not None else _v
        if _v is not None:
            self["extendiciclecolors"] = _v
        _v = arg.pop("extendpiecolors", None)
        _v = extendpiecolors if extendpiecolors is not None else _v
        if _v is not None:
            self["extendpiecolors"] = _v
        _v = arg.pop("extendsunburstcolors", None)
        _v = extendsunburstcolors if extendsunburstcolors is not None else _v
        if _v is not None:
            self["extendsunburstcolors"] = _v
        _v = arg.pop("extendtreemapcolors", None)
        _v = extendtreemapcolors if extendtreemapcolors is not None else _v
        if _v is not None:
            self["extendtreemapcolors"] = _v
        _v = arg.pop("font", None)
        _v = font if font is not None else _v
        if _v is not None:
            self["font"] = _v
        _v = arg.pop("funnelareacolorway", None)
        _v = funnelareacolorway if funnelareacolorway is not None else _v
        if _v is not None:
            self["funnelareacolorway"] = _v
        _v = arg.pop("funnelgap", None)
        _v = funnelgap if funnelgap is not None else _v
        if _v is not None:
            self["funnelgap"] = _v
        _v = arg.pop("funnelgroupgap", None)
        _v = funnelgroupgap if funnelgroupgap is not None else _v
        if _v is not None:
            self["funnelgroupgap"] = _v
        _v = arg.pop("funnelmode", None)
        _v = funnelmode if funnelmode is not None else _v
        if _v is not None:
            self["funnelmode"] = _v
        _v = arg.pop("geo", None)
        _v = geo if geo is not None else _v
        if _v is not None:
            self["geo"] = _v
        _v = arg.pop("grid", None)
        _v = grid if grid is not None else _v
        if _v is not None:
            self["grid"] = _v
        _v = arg.pop("height", None)
        _v = height if height is not None else _v
        if _v is not None:
            self["height"] = _v
        _v = arg.pop("hiddenlabels", None)
        _v = hiddenlabels if hiddenlabels is not None else _v
        if _v is not None:
            self["hiddenlabels"] = _v
        _v = arg.pop("hiddenlabelssrc", None)
        _v = hiddenlabelssrc if hiddenlabelssrc is not None else _v
        if _v is not None:
            self["hiddenlabelssrc"] = _v
        _v = arg.pop("hidesources", None)
        _v = hidesources if hidesources is not None else _v
        if _v is not None:
            self["hidesources"] = _v
        _v = arg.pop("hoverdistance", None)
        _v = hoverdistance if hoverdistance is not None else _v
        if _v is not None:
            self["hoverdistance"] = _v
        _v = arg.pop("hoverlabel", None)
        _v = hoverlabel if hoverlabel is not None else _v
        if _v is not None:
            self["hoverlabel"] = _v
        _v = arg.pop("hovermode", None)
        _v = hovermode if hovermode is not None else _v
        if _v is not None:
            self["hovermode"] = _v
        _v = arg.pop("hoversubplots", None)
        _v = hoversubplots if hoversubplots is not None else _v
        if _v is not None:
            self["hoversubplots"] = _v
        _v = arg.pop("iciclecolorway", None)
        _v = iciclecolorway if iciclecolorway is not None else _v
        if _v is not None:
            self["iciclecolorway"] = _v
        _v = arg.pop("images", None)
        _v = images if images is not None else _v
        if _v is not None:
            self["images"] = _v
        _v = arg.pop("imagedefaults", None)
        _v = imagedefaults if imagedefaults is not None else _v
        if _v is not None:
            self["imagedefaults"] = _v
        _v = arg.pop("legend", None)
        _v = legend if legend is not None else _v
        if _v is not None:
            self["legend"] = _v
        _v = arg.pop("map", None)
        _v = map if map is not None else _v
        if _v is not None:
            self["map"] = _v
        _v = arg.pop("mapbox", None)
        _v = mapbox if mapbox is not None else _v
        if _v is not None:
            self["mapbox"] = _v
        _v = arg.pop("margin", None)
        _v = margin if margin is not None else _v
        if _v is not None:
            self["margin"] = _v
        _v = arg.pop("meta", None)
        _v = meta if meta is not None else _v
        if _v is not None:
            self["meta"] = _v
        _v = arg.pop("metasrc", None)
        _v = metasrc if metasrc is not None else _v
        if _v is not None:
            self["metasrc"] = _v
        _v = arg.pop("minreducedheight", None)
        _v = minreducedheight if minreducedheight is not None else _v
        if _v is not None:
            self["minreducedheight"] = _v
        _v = arg.pop("minreducedwidth", None)
        _v = minreducedwidth if minreducedwidth is not None else _v
        if _v is not None:
            self["minreducedwidth"] = _v
        _v = arg.pop("modebar", None)
        _v = modebar if modebar is not None else _v
        if _v is not None:
            self["modebar"] = _v
        _v = arg.pop("newselection", None)
        _v = newselection if newselection is not None else _v
        if _v is not None:
            self["newselection"] = _v
        _v = arg.pop("newshape", None)
        _v = newshape if newshape is not None else _v
        if _v is not None:
            self["newshape"] = _v
        _v = arg.pop("paper_bgcolor", None)
        _v = paper_bgcolor if paper_bgcolor is not None else _v
        if _v is not None:
            self["paper_bgcolor"] = _v
        _v = arg.pop("piecolorway", None)
        _v = piecolorway if piecolorway is not None else _v
        if _v is not None:
            self["piecolorway"] = _v
        _v = arg.pop("plot_bgcolor", None)
        _v = plot_bgcolor if plot_bgcolor is not None else _v
        if _v is not None:
            self["plot_bgcolor"] = _v
        _v = arg.pop("polar", None)
        _v = polar if polar is not None else _v
        if _v is not None:
            self["polar"] = _v
        _v = arg.pop("scattergap", None)
        _v = scattergap if scattergap is not None else _v
        if _v is not None:
            self["scattergap"] = _v
        _v = arg.pop("scattermode", None)
        _v = scattermode if scattermode is not None else _v
        if _v is not None:
            self["scattermode"] = _v
        _v = arg.pop("scene", None)
        _v = scene if scene is not None else _v
        if _v is not None:
            self["scene"] = _v
        _v = arg.pop("selectdirection", None)
        _v = selectdirection if selectdirection is not None else _v
        if _v is not None:
            self["selectdirection"] = _v
        _v = arg.pop("selectionrevision", None)
        _v = selectionrevision if selectionrevision is not None else _v
        if _v is not None:
            self["selectionrevision"] = _v
        _v = arg.pop("selections", None)
        _v = selections if selections is not None else _v
        if _v is not None:
            self["selections"] = _v
        _v = arg.pop("selectiondefaults", None)
        _v = selectiondefaults if selectiondefaults is not None else _v
        if _v is not None:
            self["selectiondefaults"] = _v
        _v = arg.pop("separators", None)
        _v = separators if separators is not None else _v
        if _v is not None:
            self["separators"] = _v
        _v = arg.pop("shapes", None)
        _v = shapes if shapes is not None else _v
        if _v is not None:
            self["shapes"] = _v
        _v = arg.pop("shapedefaults", None)
        _v = shapedefaults if shapedefaults is not None else _v
        if _v is not None:
            self["shapedefaults"] = _v
        _v = arg.pop("showlegend", None)
        _v = showlegend if showlegend is not None else _v
        if _v is not None:
            self["showlegend"] = _v
        _v = arg.pop("sliders", None)
        _v = sliders if sliders is not None else _v
        if _v is not None:
            self["sliders"] = _v
        _v = arg.pop("sliderdefaults", None)
        _v = sliderdefaults if sliderdefaults is not None else _v
        if _v is not None:
            self["sliderdefaults"] = _v
        _v = arg.pop("smith", None)
        _v = smith if smith is not None else _v
        if _v is not None:
            self["smith"] = _v
        _v = arg.pop("spikedistance", None)
        _v = spikedistance if spikedistance is not None else _v
        if _v is not None:
            self["spikedistance"] = _v
        _v = arg.pop("sunburstcolorway", None)
        _v = sunburstcolorway if sunburstcolorway is not None else _v
        if _v is not None:
            self["sunburstcolorway"] = _v
        _v = arg.pop("template", None)
        _v = template if template is not None else _v
        if _v is not None:
            self["template"] = _v
        _v = arg.pop("ternary", None)
        _v = ternary if ternary is not None else _v
        if _v is not None:
            self["ternary"] = _v
        _v = arg.pop("title", None)
        _v = title if title is not None else _v
        if _v is not None:
            self["title"] = _v
        _v = arg.pop("titlefont", None)
        _v = titlefont if titlefont is not None else _v
        if _v is not None:
            self["titlefont"] = _v
        _v = arg.pop("transition", None)
        _v = transition if transition is not None else _v
        if _v is not None:
            self["transition"] = _v
        _v = arg.pop("treemapcolorway", None)
        _v = treemapcolorway if treemapcolorway is not None else _v
        if _v is not None:
            self["treemapcolorway"] = _v
        _v = arg.pop("uirevision", None)
        _v = uirevision if uirevision is not None else _v
        if _v is not None:
            self["uirevision"] = _v
        _v = arg.pop("uniformtext", None)
        _v = uniformtext if uniformtext is not None else _v
        if _v is not None:
            self["uniformtext"] = _v
        _v = arg.pop("updatemenus", None)
        _v = updatemenus if updatemenus is not None else _v
        if _v is not None:
            self["updatemenus"] = _v
        _v = arg.pop("updatemenudefaults", None)
        _v = updatemenudefaults if updatemenudefaults is not None else _v
        if _v is not None:
            self["updatemenudefaults"] = _v
        _v = arg.pop("violingap", None)
        _v = violingap if violingap is not None else _v
        if _v is not None:
            self["violingap"] = _v
        _v = arg.pop("violingroupgap", None)
        _v = violingroupgap if violingroupgap is not None else _v
        if _v is not None:
            self["violingroupgap"] = _v
        _v = arg.pop("violinmode", None)
        _v = violinmode if violinmode is not None else _v
        if _v is not None:
            self["violinmode"] = _v
        _v = arg.pop("waterfallgap", None)
        _v = waterfallgap if waterfallgap is not None else _v
        if _v is not None:
            self["waterfallgap"] = _v
        _v = arg.pop("waterfallgroupgap", None)
        _v = waterfallgroupgap if waterfallgroupgap is not None else _v
        if _v is not None:
            self["waterfallgroupgap"] = _v
        _v = arg.pop("waterfallmode", None)
        _v = waterfallmode if waterfallmode is not None else _v
        if _v is not None:
            self["waterfallmode"] = _v
        _v = arg.pop("width", None)
        _v = width if width is not None else _v
        if _v is not None:
            self["width"] = _v
        _v = arg.pop("xaxis", None)
        _v = xaxis if xaxis is not None else _v
        if _v is not None:
            self["xaxis"] = _v
        _v = arg.pop("yaxis", None)
        _v = yaxis if yaxis is not None else _v
        if _v is not None:
            self["yaxis"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
