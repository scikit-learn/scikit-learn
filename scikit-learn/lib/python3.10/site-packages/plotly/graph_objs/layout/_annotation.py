#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Annotation(_BaseLayoutHierarchyType):
    _parent_path_str = "layout"
    _path_str = "layout.annotation"
    _valid_props = {
        "align",
        "arrowcolor",
        "arrowhead",
        "arrowside",
        "arrowsize",
        "arrowwidth",
        "ax",
        "axref",
        "ay",
        "ayref",
        "bgcolor",
        "bordercolor",
        "borderpad",
        "borderwidth",
        "captureevents",
        "clicktoshow",
        "font",
        "height",
        "hoverlabel",
        "hovertext",
        "name",
        "opacity",
        "showarrow",
        "standoff",
        "startarrowhead",
        "startarrowsize",
        "startstandoff",
        "templateitemname",
        "text",
        "textangle",
        "valign",
        "visible",
        "width",
        "x",
        "xanchor",
        "xclick",
        "xref",
        "xshift",
        "y",
        "yanchor",
        "yclick",
        "yref",
        "yshift",
    }

    @property
    def align(self):
        """
        Sets the horizontal alignment of the `text` within the box. Has
        an effect only if `text` spans two or more lines (i.e. `text`
        contains one or more <br> HTML tags) or if an explicit width is
        set to override the text width.

        The 'align' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['left', 'center', 'right']

        Returns
        -------
        Any
        """
        return self["align"]

    @align.setter
    def align(self, val):
        self["align"] = val

    @property
    def arrowcolor(self):
        """
        Sets the color of the annotation arrow.

        The 'arrowcolor' property is a color and may be specified as:
          - A hex string (e.g. '#ff0000')
          - An rgb/rgba string (e.g. 'rgb(255,0,0)')
          - An hsl/hsla string (e.g. 'hsl(0,100%,50%)')
          - An hsv/hsva string (e.g. 'hsv(0,100%,100%)')
          - A named CSS color: see https://plotly.com/python/css-colors/ for a list

        Returns
        -------
        str
        """
        return self["arrowcolor"]

    @arrowcolor.setter
    def arrowcolor(self, val):
        self["arrowcolor"] = val

    @property
    def arrowhead(self):
        """
        Sets the end annotation arrow head style.

        The 'arrowhead' property is a integer and may be specified as:
          - An int (or float that will be cast to an int)
            in the interval [0, 8]

        Returns
        -------
        int
        """
        return self["arrowhead"]

    @arrowhead.setter
    def arrowhead(self, val):
        self["arrowhead"] = val

    @property
    def arrowside(self):
        """
        Sets the annotation arrow head position.

        The 'arrowside' property is a flaglist and may be specified
        as a string containing:
          - Any combination of ['end', 'start'] joined with '+' characters
            (e.g. 'end+start')
            OR exactly one of ['none'] (e.g. 'none')

        Returns
        -------
        Any
        """
        return self["arrowside"]

    @arrowside.setter
    def arrowside(self, val):
        self["arrowside"] = val

    @property
    def arrowsize(self):
        """
        Sets the size of the end annotation arrow head, relative to
        `arrowwidth`. A value of 1 (default) gives a head about 3x as
        wide as the line.

        The 'arrowsize' property is a number and may be specified as:
          - An int or float in the interval [0.3, inf]

        Returns
        -------
        int|float
        """
        return self["arrowsize"]

    @arrowsize.setter
    def arrowsize(self, val):
        self["arrowsize"] = val

    @property
    def arrowwidth(self):
        """
        Sets the width (in px) of annotation arrow line.

        The 'arrowwidth' property is a number and may be specified as:
          - An int or float in the interval [0.1, inf]

        Returns
        -------
        int|float
        """
        return self["arrowwidth"]

    @arrowwidth.setter
    def arrowwidth(self, val):
        self["arrowwidth"] = val

    @property
    def ax(self):
        """
        Sets the x component of the arrow tail about the arrow head. If
        `axref` is `pixel`, a positive (negative) component corresponds
        to an arrow pointing from right to left (left to right). If
        `axref` is not `pixel` and is exactly the same as `xref`, this
        is an absolute value on that axis, like `x`, specified in the
        same coordinates as `xref`.

        The 'ax' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["ax"]

    @ax.setter
    def ax(self, val):
        self["ax"] = val

    @property
    def axref(self):
        """
        Indicates in what coordinates the tail of the annotation
        (ax,ay) is specified. If set to a x axis id (e.g. "x" or "x2"),
        the `x` position refers to a x coordinate. If set to "paper",
        the `x` position refers to the distance from the left of the
        plotting area in normalized coordinates where 0 (1) corresponds
        to the left (right). If set to a x axis ID followed by "domain"
        (separated by a space), the position behaves like for "paper",
        but refers to the distance in fractions of the domain length
        from the left of the domain of that axis: e.g., *x2 domain*
        refers to the domain of the second x  axis and a x position of
        0.5 refers to the point between the left and the right of the
        domain of the second x axis. In order for absolute positioning
        of the arrow to work, "axref" must be exactly the same as
        "xref", otherwise "axref" will revert to "pixel" (explained
        next). For relative positioning, "axref" can be set to "pixel",
        in which case the "ax" value is specified in pixels relative to
        "x". Absolute positioning is useful for trendline annotations
        which should continue to indicate the correct trend when
        zoomed. Relative positioning is useful for specifying the text
        offset for an annotated point.

        The 'axref' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['pixel']
          - A string that matches one of the following regular expressions:
                ['^x([2-9]|[1-9][0-9]+)?( domain)?$']

        Returns
        -------
        Any
        """
        return self["axref"]

    @axref.setter
    def axref(self, val):
        self["axref"] = val

    @property
    def ay(self):
        """
        Sets the y component of the arrow tail about the arrow head. If
        `ayref` is `pixel`, a positive (negative) component corresponds
        to an arrow pointing from bottom to top (top to bottom). If
        `ayref` is not `pixel` and is exactly the same as `yref`, this
        is an absolute value on that axis, like `y`, specified in the
        same coordinates as `yref`.

        The 'ay' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["ay"]

    @ay.setter
    def ay(self, val):
        self["ay"] = val

    @property
    def ayref(self):
        """
        Indicates in what coordinates the tail of the annotation
        (ax,ay) is specified. If set to a y axis id (e.g. "y" or "y2"),
        the `y` position refers to a y coordinate. If set to "paper",
        the `y` position refers to the distance from the bottom of the
        plotting area in normalized coordinates where 0 (1) corresponds
        to the bottom (top). If set to a y axis ID followed by "domain"
        (separated by a space), the position behaves like for "paper",
        but refers to the distance in fractions of the domain length
        from the bottom of the domain of that axis: e.g., *y2 domain*
        refers to the domain of the second y  axis and a y position of
        0.5 refers to the point between the bottom and the top of the
        domain of the second y axis. In order for absolute positioning
        of the arrow to work, "ayref" must be exactly the same as
        "yref", otherwise "ayref" will revert to "pixel" (explained
        next). For relative positioning, "ayref" can be set to "pixel",
        in which case the "ay" value is specified in pixels relative to
        "y". Absolute positioning is useful for trendline annotations
        which should continue to indicate the correct trend when
        zoomed. Relative positioning is useful for specifying the text
        offset for an annotated point.

        The 'ayref' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['pixel']
          - A string that matches one of the following regular expressions:
                ['^y([2-9]|[1-9][0-9]+)?( domain)?$']

        Returns
        -------
        Any
        """
        return self["ayref"]

    @ayref.setter
    def ayref(self, val):
        self["ayref"] = val

    @property
    def bgcolor(self):
        """
        Sets the background color of the annotation.

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
        Sets the color of the border enclosing the annotation `text`.

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
    def borderpad(self):
        """
        Sets the padding (in px) between the `text` and the enclosing
        border.

        The 'borderpad' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["borderpad"]

    @borderpad.setter
    def borderpad(self, val):
        self["borderpad"] = val

    @property
    def borderwidth(self):
        """
        Sets the width (in px) of the border enclosing the annotation
        `text`.

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
    def captureevents(self):
        """
        Determines whether the annotation text box captures mouse move
        and click events, or allows those events to pass through to
        data points in the plot that may be behind the annotation. By
        default `captureevents` is False unless `hovertext` is
        provided. If you use the event `plotly_clickannotation` without
        `hovertext` you must explicitly enable `captureevents`.

        The 'captureevents' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["captureevents"]

    @captureevents.setter
    def captureevents(self, val):
        self["captureevents"] = val

    @property
    def clicktoshow(self):
        """
        Makes this annotation respond to clicks on the plot. If you
        click a data point that exactly matches the `x` and `y` values
        of this annotation, and it is hidden (visible: false), it will
        appear. In "onoff" mode, you must click the same point again to
        make it disappear, so if you click multiple points, you can
        show multiple annotations. In "onout" mode, a click anywhere
        else in the plot (on another data point or not) will hide this
        annotation. If you need to show/hide this annotation in
        response to different `x` or `y` values, you can set `xclick`
        and/or `yclick`. This is useful for example to label the side
        of a bar. To label markers though, `standoff` is preferred over
        `xclick` and `yclick`.

        The 'clicktoshow' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                [False, 'onoff', 'onout']

        Returns
        -------
        Any
        """
        return self["clicktoshow"]

    @clicktoshow.setter
    def clicktoshow(self, val):
        self["clicktoshow"] = val

    @property
    def font(self):
        """
        Sets the annotation text font.

        The 'font' property is an instance of Font
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.annotation.Font`
          - A dict of string/value properties that will be passed
            to the Font constructor

        Returns
        -------
        plotly.graph_objs.layout.annotation.Font
        """
        return self["font"]

    @font.setter
    def font(self, val):
        self["font"] = val

    @property
    def height(self):
        """
        Sets an explicit height for the text box. null (default) lets
        the text set the box height. Taller text will be clipped.

        The 'height' property is a number and may be specified as:
          - An int or float in the interval [1, inf]

        Returns
        -------
        int|float
        """
        return self["height"]

    @height.setter
    def height(self, val):
        self["height"] = val

    @property
    def hoverlabel(self):
        """
        The 'hoverlabel' property is an instance of Hoverlabel
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.annotation.Hoverlabel`
          - A dict of string/value properties that will be passed
            to the Hoverlabel constructor

        Returns
        -------
        plotly.graph_objs.layout.annotation.Hoverlabel
        """
        return self["hoverlabel"]

    @hoverlabel.setter
    def hoverlabel(self, val):
        self["hoverlabel"] = val

    @property
    def hovertext(self):
        """
        Sets text to appear when hovering over this annotation. If
        omitted or blank, no hover label will appear.

        The 'hovertext' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["hovertext"]

    @hovertext.setter
    def hovertext(self, val):
        self["hovertext"] = val

    @property
    def name(self):
        """
        When used in a template, named items are created in the output
        figure in addition to any items the figure already has in this
        array. You can modify these items in the output figure by
        making your own item with `templateitemname` matching this
        `name` alongside your modifications (including `visible: false`
        or `enabled: false` to hide it). Has no effect outside of a
        template.

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

    @property
    def opacity(self):
        """
        Sets the opacity of the annotation (text + arrow).

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

    @property
    def showarrow(self):
        """
        Determines whether or not the annotation is drawn with an
        arrow. If True, `text` is placed near the arrow's tail. If
        False, `text` lines up with the `x` and `y` provided.

        The 'showarrow' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["showarrow"]

    @showarrow.setter
    def showarrow(self, val):
        self["showarrow"] = val

    @property
    def standoff(self):
        """
        Sets a distance, in pixels, to move the end arrowhead away from
        the position it is pointing at, for example to point at the
        edge of a marker independent of zoom. Note that this shortens
        the arrow from the `ax` / `ay` vector, in contrast to `xshift`
        / `yshift` which moves everything by this amount.

        The 'standoff' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["standoff"]

    @standoff.setter
    def standoff(self, val):
        self["standoff"] = val

    @property
    def startarrowhead(self):
        """
        Sets the start annotation arrow head style.

        The 'startarrowhead' property is a integer and may be specified as:
          - An int (or float that will be cast to an int)
            in the interval [0, 8]

        Returns
        -------
        int
        """
        return self["startarrowhead"]

    @startarrowhead.setter
    def startarrowhead(self, val):
        self["startarrowhead"] = val

    @property
    def startarrowsize(self):
        """
        Sets the size of the start annotation arrow head, relative to
        `arrowwidth`. A value of 1 (default) gives a head about 3x as
        wide as the line.

        The 'startarrowsize' property is a number and may be specified as:
          - An int or float in the interval [0.3, inf]

        Returns
        -------
        int|float
        """
        return self["startarrowsize"]

    @startarrowsize.setter
    def startarrowsize(self, val):
        self["startarrowsize"] = val

    @property
    def startstandoff(self):
        """
        Sets a distance, in pixels, to move the start arrowhead away
        from the position it is pointing at, for example to point at
        the edge of a marker independent of zoom. Note that this
        shortens the arrow from the `ax` / `ay` vector, in contrast to
        `xshift` / `yshift` which moves everything by this amount.

        The 'startstandoff' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["startstandoff"]

    @startstandoff.setter
    def startstandoff(self, val):
        self["startstandoff"] = val

    @property
    def templateitemname(self):
        """
        Used to refer to a named item in this array in the template.
        Named items from the template will be created even without a
        matching item in the input figure, but you can modify one by
        making an item with `templateitemname` matching its `name`,
        alongside your modifications (including `visible: false` or
        `enabled: false` to hide it). If there is no template or no
        matching item, this item will be hidden unless you explicitly
        show it with `visible: true`.

        The 'templateitemname' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["templateitemname"]

    @templateitemname.setter
    def templateitemname(self, val):
        self["templateitemname"] = val

    @property
    def text(self):
        """
        Sets the text associated with this annotation. Plotly uses a
        subset of HTML tags to do things like newline (`<br>`), bold
        (`<b></b>`), italics (`<i></i>`), hyperlinks (`<a
        href='...'></a>`). Tags `<em>`, `<sup>`, `<sub>`, `<s>`, `<u>`,
        and `<span>` are also supported.

        The 'text' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["text"]

    @text.setter
    def text(self, val):
        self["text"] = val

    @property
    def textangle(self):
        """
        Sets the angle at which the `text` is drawn with respect to the
        horizontal.

        The 'textangle' property is a angle (in degrees) that may be
        specified as a number between -180 and 180.
        Numeric values outside this range are converted to the equivalent value
        (e.g. 270 is converted to -90).

        Returns
        -------
        int|float
        """
        return self["textangle"]

    @textangle.setter
    def textangle(self, val):
        self["textangle"] = val

    @property
    def valign(self):
        """
        Sets the vertical alignment of the `text` within the box. Has
        an effect only if an explicit height is set to override the
        text height.

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
        Determines whether or not this annotation is visible.

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
    def width(self):
        """
        Sets an explicit width for the text box. null (default) lets
        the text set the box width. Wider text will be clipped. There
        is no automatic wrapping; use <br> to start a new line.

        The 'width' property is a number and may be specified as:
          - An int or float in the interval [1, inf]

        Returns
        -------
        int|float
        """
        return self["width"]

    @width.setter
    def width(self, val):
        self["width"] = val

    @property
    def x(self):
        """
        Sets the annotation's x position. If the axis `type` is "log",
        then you must take the log of your desired range. If the axis
        `type` is "date", it should be date strings, like date data,
        though Date objects and unix milliseconds will be accepted and
        converted to strings. If the axis `type` is "category", it
        should be numbers, using the scale where each category is
        assigned a serial number from zero in the order it appears.

        The 'x' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["x"]

    @x.setter
    def x(self, val):
        self["x"] = val

    @property
    def xanchor(self):
        """
        Sets the text box's horizontal position anchor This anchor
        binds the `x` position to the "left", "center" or "right" of
        the annotation. For example, if `x` is set to 1, `xref` to
        "paper" and `xanchor` to "right" then the right-most portion of
        the annotation lines up with the right-most edge of the
        plotting area. If "auto", the anchor is equivalent to "center"
        for data-referenced annotations or if there is an arrow,
        whereas for paper-referenced with no arrow, the anchor picked
        corresponds to the closest side.

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
    def xclick(self):
        """
        Toggle this annotation when clicking a data point whose `x`
        value is `xclick` rather than the annotation's `x` value.

        The 'xclick' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["xclick"]

    @xclick.setter
    def xclick(self, val):
        self["xclick"] = val

    @property
    def xref(self):
        """
        Sets the annotation's x coordinate axis. If set to a x axis id
        (e.g. "x" or "x2"), the `x` position refers to a x coordinate.
        If set to "paper", the `x` position refers to the distance from
        the left of the plotting area in normalized coordinates where 0
        (1) corresponds to the left (right). If set to a x axis ID
        followed by "domain" (separated by a space), the position
        behaves like for "paper", but refers to the distance in
        fractions of the domain length from the left of the domain of
        that axis: e.g., *x2 domain* refers to the domain of the second
        x  axis and a x position of 0.5 refers to the point between the
        left and the right of the domain of the second x axis.

        The 'xref' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['paper']
          - A string that matches one of the following regular expressions:
                ['^x([2-9]|[1-9][0-9]+)?( domain)?$']

        Returns
        -------
        Any
        """
        return self["xref"]

    @xref.setter
    def xref(self, val):
        self["xref"] = val

    @property
    def xshift(self):
        """
        Shifts the position of the whole annotation and arrow to the
        right (positive) or left (negative) by this many pixels.

        The 'xshift' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["xshift"]

    @xshift.setter
    def xshift(self, val):
        self["xshift"] = val

    @property
    def y(self):
        """
        Sets the annotation's y position. If the axis `type` is "log",
        then you must take the log of your desired range. If the axis
        `type` is "date", it should be date strings, like date data,
        though Date objects and unix milliseconds will be accepted and
        converted to strings. If the axis `type` is "category", it
        should be numbers, using the scale where each category is
        assigned a serial number from zero in the order it appears.

        The 'y' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["y"]

    @y.setter
    def y(self, val):
        self["y"] = val

    @property
    def yanchor(self):
        """
        Sets the text box's vertical position anchor This anchor binds
        the `y` position to the "top", "middle" or "bottom" of the
        annotation. For example, if `y` is set to 1, `yref` to "paper"
        and `yanchor` to "top" then the top-most portion of the
        annotation lines up with the top-most edge of the plotting
        area. If "auto", the anchor is equivalent to "middle" for data-
        referenced annotations or if there is an arrow, whereas for
        paper-referenced with no arrow, the anchor picked corresponds
        to the closest side.

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
    def yclick(self):
        """
        Toggle this annotation when clicking a data point whose `y`
        value is `yclick` rather than the annotation's `y` value.

        The 'yclick' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["yclick"]

    @yclick.setter
    def yclick(self, val):
        self["yclick"] = val

    @property
    def yref(self):
        """
        Sets the annotation's y coordinate axis. If set to a y axis id
        (e.g. "y" or "y2"), the `y` position refers to a y coordinate.
        If set to "paper", the `y` position refers to the distance from
        the bottom of the plotting area in normalized coordinates where
        0 (1) corresponds to the bottom (top). If set to a y axis ID
        followed by "domain" (separated by a space), the position
        behaves like for "paper", but refers to the distance in
        fractions of the domain length from the bottom of the domain of
        that axis: e.g., *y2 domain* refers to the domain of the second
        y  axis and a y position of 0.5 refers to the point between the
        bottom and the top of the domain of the second y axis.

        The 'yref' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['paper']
          - A string that matches one of the following regular expressions:
                ['^y([2-9]|[1-9][0-9]+)?( domain)?$']

        Returns
        -------
        Any
        """
        return self["yref"]

    @yref.setter
    def yref(self, val):
        self["yref"] = val

    @property
    def yshift(self):
        """
        Shifts the position of the whole annotation and arrow up
        (positive) or down (negative) by this many pixels.

        The 'yshift' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["yshift"]

    @yshift.setter
    def yshift(self, val):
        self["yshift"] = val

    @property
    def _prop_descriptions(self):
        return """\
        align
            Sets the horizontal alignment of the `text` within the
            box. Has an effect only if `text` spans two or more
            lines (i.e. `text` contains one or more <br> HTML tags)
            or if an explicit width is set to override the text
            width.
        arrowcolor
            Sets the color of the annotation arrow.
        arrowhead
            Sets the end annotation arrow head style.
        arrowside
            Sets the annotation arrow head position.
        arrowsize
            Sets the size of the end annotation arrow head,
            relative to `arrowwidth`. A value of 1 (default) gives
            a head about 3x as wide as the line.
        arrowwidth
            Sets the width (in px) of annotation arrow line.
        ax
            Sets the x component of the arrow tail about the arrow
            head. If `axref` is `pixel`, a positive (negative)
            component corresponds to an arrow pointing from right
            to left (left to right). If `axref` is not `pixel` and
            is exactly the same as `xref`, this is an absolute
            value on that axis, like `x`, specified in the same
            coordinates as `xref`.
        axref
            Indicates in what coordinates the tail of the
            annotation (ax,ay) is specified. If set to a x axis id
            (e.g. "x" or "x2"), the `x` position refers to a x
            coordinate. If set to "paper", the `x` position refers
            to the distance from the left of the plotting area in
            normalized coordinates where 0 (1) corresponds to the
            left (right). If set to a x axis ID followed by
            "domain" (separated by a space), the position behaves
            like for "paper", but refers to the distance in
            fractions of the domain length from the left of the
            domain of that axis: e.g., *x2 domain* refers to the
            domain of the second x  axis and a x position of 0.5
            refers to the point between the left and the right of
            the domain of the second x axis. In order for absolute
            positioning of the arrow to work, "axref" must be
            exactly the same as "xref", otherwise "axref" will
            revert to "pixel" (explained next). For relative
            positioning, "axref" can be set to "pixel", in which
            case the "ax" value is specified in pixels relative to
            "x". Absolute positioning is useful for trendline
            annotations which should continue to indicate the
            correct trend when zoomed. Relative positioning is
            useful for specifying the text offset for an annotated
            point.
        ay
            Sets the y component of the arrow tail about the arrow
            head. If `ayref` is `pixel`, a positive (negative)
            component corresponds to an arrow pointing from bottom
            to top (top to bottom). If `ayref` is not `pixel` and
            is exactly the same as `yref`, this is an absolute
            value on that axis, like `y`, specified in the same
            coordinates as `yref`.
        ayref
            Indicates in what coordinates the tail of the
            annotation (ax,ay) is specified. If set to a y axis id
            (e.g. "y" or "y2"), the `y` position refers to a y
            coordinate. If set to "paper", the `y` position refers
            to the distance from the bottom of the plotting area in
            normalized coordinates where 0 (1) corresponds to the
            bottom (top). If set to a y axis ID followed by
            "domain" (separated by a space), the position behaves
            like for "paper", but refers to the distance in
            fractions of the domain length from the bottom of the
            domain of that axis: e.g., *y2 domain* refers to the
            domain of the second y  axis and a y position of 0.5
            refers to the point between the bottom and the top of
            the domain of the second y axis. In order for absolute
            positioning of the arrow to work, "ayref" must be
            exactly the same as "yref", otherwise "ayref" will
            revert to "pixel" (explained next). For relative
            positioning, "ayref" can be set to "pixel", in which
            case the "ay" value is specified in pixels relative to
            "y". Absolute positioning is useful for trendline
            annotations which should continue to indicate the
            correct trend when zoomed. Relative positioning is
            useful for specifying the text offset for an annotated
            point.
        bgcolor
            Sets the background color of the annotation.
        bordercolor
            Sets the color of the border enclosing the annotation
            `text`.
        borderpad
            Sets the padding (in px) between the `text` and the
            enclosing border.
        borderwidth
            Sets the width (in px) of the border enclosing the
            annotation `text`.
        captureevents
            Determines whether the annotation text box captures
            mouse move and click events, or allows those events to
            pass through to data points in the plot that may be
            behind the annotation. By default `captureevents` is
            False unless `hovertext` is provided. If you use the
            event `plotly_clickannotation` without `hovertext` you
            must explicitly enable `captureevents`.
        clicktoshow
            Makes this annotation respond to clicks on the plot. If
            you click a data point that exactly matches the `x` and
            `y` values of this annotation, and it is hidden
            (visible: false), it will appear. In "onoff" mode, you
            must click the same point again to make it disappear,
            so if you click multiple points, you can show multiple
            annotations. In "onout" mode, a click anywhere else in
            the plot (on another data point or not) will hide this
            annotation. If you need to show/hide this annotation in
            response to different `x` or `y` values, you can set
            `xclick` and/or `yclick`. This is useful for example to
            label the side of a bar. To label markers though,
            `standoff` is preferred over `xclick` and `yclick`.
        font
            Sets the annotation text font.
        height
            Sets an explicit height for the text box. null
            (default) lets the text set the box height. Taller text
            will be clipped.
        hoverlabel
            :class:`plotly.graph_objects.layout.annotation.Hoverlab
            el` instance or dict with compatible properties
        hovertext
            Sets text to appear when hovering over this annotation.
            If omitted or blank, no hover label will appear.
        name
            When used in a template, named items are created in the
            output figure in addition to any items the figure
            already has in this array. You can modify these items
            in the output figure by making your own item with
            `templateitemname` matching this `name` alongside your
            modifications (including `visible: false` or `enabled:
            false` to hide it). Has no effect outside of a
            template.
        opacity
            Sets the opacity of the annotation (text + arrow).
        showarrow
            Determines whether or not the annotation is drawn with
            an arrow. If True, `text` is placed near the arrow's
            tail. If False, `text` lines up with the `x` and `y`
            provided.
        standoff
            Sets a distance, in pixels, to move the end arrowhead
            away from the position it is pointing at, for example
            to point at the edge of a marker independent of zoom.
            Note that this shortens the arrow from the `ax` / `ay`
            vector, in contrast to `xshift` / `yshift` which moves
            everything by this amount.
        startarrowhead
            Sets the start annotation arrow head style.
        startarrowsize
            Sets the size of the start annotation arrow head,
            relative to `arrowwidth`. A value of 1 (default) gives
            a head about 3x as wide as the line.
        startstandoff
            Sets a distance, in pixels, to move the start arrowhead
            away from the position it is pointing at, for example
            to point at the edge of a marker independent of zoom.
            Note that this shortens the arrow from the `ax` / `ay`
            vector, in contrast to `xshift` / `yshift` which moves
            everything by this amount.
        templateitemname
            Used to refer to a named item in this array in the
            template. Named items from the template will be created
            even without a matching item in the input figure, but
            you can modify one by making an item with
            `templateitemname` matching its `name`, alongside your
            modifications (including `visible: false` or `enabled:
            false` to hide it). If there is no template or no
            matching item, this item will be hidden unless you
            explicitly show it with `visible: true`.
        text
            Sets the text associated with this annotation. Plotly
            uses a subset of HTML tags to do things like newline
            (`<br>`), bold (`<b></b>`), italics (`<i></i>`),
            hyperlinks (`<a href='...'></a>`). Tags `<em>`,
            `<sup>`, `<sub>`, `<s>`, `<u>`, and `<span>` are also
            supported.
        textangle
            Sets the angle at which the `text` is drawn with
            respect to the horizontal.
        valign
            Sets the vertical alignment of the `text` within the
            box. Has an effect only if an explicit height is set to
            override the text height.
        visible
            Determines whether or not this annotation is visible.
        width
            Sets an explicit width for the text box. null (default)
            lets the text set the box width. Wider text will be
            clipped. There is no automatic wrapping; use <br> to
            start a new line.
        x
            Sets the annotation's x position. If the axis `type` is
            "log", then you must take the log of your desired
            range. If the axis `type` is "date", it should be date
            strings, like date data, though Date objects and unix
            milliseconds will be accepted and converted to strings.
            If the axis `type` is "category", it should be numbers,
            using the scale where each category is assigned a
            serial number from zero in the order it appears.
        xanchor
            Sets the text box's horizontal position anchor This
            anchor binds the `x` position to the "left", "center"
            or "right" of the annotation. For example, if `x` is
            set to 1, `xref` to "paper" and `xanchor` to "right"
            then the right-most portion of the annotation lines up
            with the right-most edge of the plotting area. If
            "auto", the anchor is equivalent to "center" for data-
            referenced annotations or if there is an arrow, whereas
            for paper-referenced with no arrow, the anchor picked
            corresponds to the closest side.
        xclick
            Toggle this annotation when clicking a data point whose
            `x` value is `xclick` rather than the annotation's `x`
            value.
        xref
            Sets the annotation's x coordinate axis. If set to a x
            axis id (e.g. "x" or "x2"), the `x` position refers to
            a x coordinate. If set to "paper", the `x` position
            refers to the distance from the left of the plotting
            area in normalized coordinates where 0 (1) corresponds
            to the left (right). If set to a x axis ID followed by
            "domain" (separated by a space), the position behaves
            like for "paper", but refers to the distance in
            fractions of the domain length from the left of the
            domain of that axis: e.g., *x2 domain* refers to the
            domain of the second x  axis and a x position of 0.5
            refers to the point between the left and the right of
            the domain of the second x axis.
        xshift
            Shifts the position of the whole annotation and arrow
            to the right (positive) or left (negative) by this many
            pixels.
        y
            Sets the annotation's y position. If the axis `type` is
            "log", then you must take the log of your desired
            range. If the axis `type` is "date", it should be date
            strings, like date data, though Date objects and unix
            milliseconds will be accepted and converted to strings.
            If the axis `type` is "category", it should be numbers,
            using the scale where each category is assigned a
            serial number from zero in the order it appears.
        yanchor
            Sets the text box's vertical position anchor This
            anchor binds the `y` position to the "top", "middle" or
            "bottom" of the annotation. For example, if `y` is set
            to 1, `yref` to "paper" and `yanchor` to "top" then the
            top-most portion of the annotation lines up with the
            top-most edge of the plotting area. If "auto", the
            anchor is equivalent to "middle" for data-referenced
            annotations or if there is an arrow, whereas for paper-
            referenced with no arrow, the anchor picked corresponds
            to the closest side.
        yclick
            Toggle this annotation when clicking a data point whose
            `y` value is `yclick` rather than the annotation's `y`
            value.
        yref
            Sets the annotation's y coordinate axis. If set to a y
            axis id (e.g. "y" or "y2"), the `y` position refers to
            a y coordinate. If set to "paper", the `y` position
            refers to the distance from the bottom of the plotting
            area in normalized coordinates where 0 (1) corresponds
            to the bottom (top). If set to a y axis ID followed by
            "domain" (separated by a space), the position behaves
            like for "paper", but refers to the distance in
            fractions of the domain length from the bottom of the
            domain of that axis: e.g., *y2 domain* refers to the
            domain of the second y  axis and a y position of 0.5
            refers to the point between the bottom and the top of
            the domain of the second y axis.
        yshift
            Shifts the position of the whole annotation and arrow
            up (positive) or down (negative) by this many pixels.
        """

    def __init__(
        self,
        arg=None,
        align=None,
        arrowcolor=None,
        arrowhead=None,
        arrowside=None,
        arrowsize=None,
        arrowwidth=None,
        ax=None,
        axref=None,
        ay=None,
        ayref=None,
        bgcolor=None,
        bordercolor=None,
        borderpad=None,
        borderwidth=None,
        captureevents=None,
        clicktoshow=None,
        font=None,
        height=None,
        hoverlabel=None,
        hovertext=None,
        name=None,
        opacity=None,
        showarrow=None,
        standoff=None,
        startarrowhead=None,
        startarrowsize=None,
        startstandoff=None,
        templateitemname=None,
        text=None,
        textangle=None,
        valign=None,
        visible=None,
        width=None,
        x=None,
        xanchor=None,
        xclick=None,
        xref=None,
        xshift=None,
        y=None,
        yanchor=None,
        yclick=None,
        yref=None,
        yshift=None,
        **kwargs,
    ):
        """
        Construct a new Annotation object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.layout.Annotation`
        align
            Sets the horizontal alignment of the `text` within the
            box. Has an effect only if `text` spans two or more
            lines (i.e. `text` contains one or more <br> HTML tags)
            or if an explicit width is set to override the text
            width.
        arrowcolor
            Sets the color of the annotation arrow.
        arrowhead
            Sets the end annotation arrow head style.
        arrowside
            Sets the annotation arrow head position.
        arrowsize
            Sets the size of the end annotation arrow head,
            relative to `arrowwidth`. A value of 1 (default) gives
            a head about 3x as wide as the line.
        arrowwidth
            Sets the width (in px) of annotation arrow line.
        ax
            Sets the x component of the arrow tail about the arrow
            head. If `axref` is `pixel`, a positive (negative)
            component corresponds to an arrow pointing from right
            to left (left to right). If `axref` is not `pixel` and
            is exactly the same as `xref`, this is an absolute
            value on that axis, like `x`, specified in the same
            coordinates as `xref`.
        axref
            Indicates in what coordinates the tail of the
            annotation (ax,ay) is specified. If set to a x axis id
            (e.g. "x" or "x2"), the `x` position refers to a x
            coordinate. If set to "paper", the `x` position refers
            to the distance from the left of the plotting area in
            normalized coordinates where 0 (1) corresponds to the
            left (right). If set to a x axis ID followed by
            "domain" (separated by a space), the position behaves
            like for "paper", but refers to the distance in
            fractions of the domain length from the left of the
            domain of that axis: e.g., *x2 domain* refers to the
            domain of the second x  axis and a x position of 0.5
            refers to the point between the left and the right of
            the domain of the second x axis. In order for absolute
            positioning of the arrow to work, "axref" must be
            exactly the same as "xref", otherwise "axref" will
            revert to "pixel" (explained next). For relative
            positioning, "axref" can be set to "pixel", in which
            case the "ax" value is specified in pixels relative to
            "x". Absolute positioning is useful for trendline
            annotations which should continue to indicate the
            correct trend when zoomed. Relative positioning is
            useful for specifying the text offset for an annotated
            point.
        ay
            Sets the y component of the arrow tail about the arrow
            head. If `ayref` is `pixel`, a positive (negative)
            component corresponds to an arrow pointing from bottom
            to top (top to bottom). If `ayref` is not `pixel` and
            is exactly the same as `yref`, this is an absolute
            value on that axis, like `y`, specified in the same
            coordinates as `yref`.
        ayref
            Indicates in what coordinates the tail of the
            annotation (ax,ay) is specified. If set to a y axis id
            (e.g. "y" or "y2"), the `y` position refers to a y
            coordinate. If set to "paper", the `y` position refers
            to the distance from the bottom of the plotting area in
            normalized coordinates where 0 (1) corresponds to the
            bottom (top). If set to a y axis ID followed by
            "domain" (separated by a space), the position behaves
            like for "paper", but refers to the distance in
            fractions of the domain length from the bottom of the
            domain of that axis: e.g., *y2 domain* refers to the
            domain of the second y  axis and a y position of 0.5
            refers to the point between the bottom and the top of
            the domain of the second y axis. In order for absolute
            positioning of the arrow to work, "ayref" must be
            exactly the same as "yref", otherwise "ayref" will
            revert to "pixel" (explained next). For relative
            positioning, "ayref" can be set to "pixel", in which
            case the "ay" value is specified in pixels relative to
            "y". Absolute positioning is useful for trendline
            annotations which should continue to indicate the
            correct trend when zoomed. Relative positioning is
            useful for specifying the text offset for an annotated
            point.
        bgcolor
            Sets the background color of the annotation.
        bordercolor
            Sets the color of the border enclosing the annotation
            `text`.
        borderpad
            Sets the padding (in px) between the `text` and the
            enclosing border.
        borderwidth
            Sets the width (in px) of the border enclosing the
            annotation `text`.
        captureevents
            Determines whether the annotation text box captures
            mouse move and click events, or allows those events to
            pass through to data points in the plot that may be
            behind the annotation. By default `captureevents` is
            False unless `hovertext` is provided. If you use the
            event `plotly_clickannotation` without `hovertext` you
            must explicitly enable `captureevents`.
        clicktoshow
            Makes this annotation respond to clicks on the plot. If
            you click a data point that exactly matches the `x` and
            `y` values of this annotation, and it is hidden
            (visible: false), it will appear. In "onoff" mode, you
            must click the same point again to make it disappear,
            so if you click multiple points, you can show multiple
            annotations. In "onout" mode, a click anywhere else in
            the plot (on another data point or not) will hide this
            annotation. If you need to show/hide this annotation in
            response to different `x` or `y` values, you can set
            `xclick` and/or `yclick`. This is useful for example to
            label the side of a bar. To label markers though,
            `standoff` is preferred over `xclick` and `yclick`.
        font
            Sets the annotation text font.
        height
            Sets an explicit height for the text box. null
            (default) lets the text set the box height. Taller text
            will be clipped.
        hoverlabel
            :class:`plotly.graph_objects.layout.annotation.Hoverlab
            el` instance or dict with compatible properties
        hovertext
            Sets text to appear when hovering over this annotation.
            If omitted or blank, no hover label will appear.
        name
            When used in a template, named items are created in the
            output figure in addition to any items the figure
            already has in this array. You can modify these items
            in the output figure by making your own item with
            `templateitemname` matching this `name` alongside your
            modifications (including `visible: false` or `enabled:
            false` to hide it). Has no effect outside of a
            template.
        opacity
            Sets the opacity of the annotation (text + arrow).
        showarrow
            Determines whether or not the annotation is drawn with
            an arrow. If True, `text` is placed near the arrow's
            tail. If False, `text` lines up with the `x` and `y`
            provided.
        standoff
            Sets a distance, in pixels, to move the end arrowhead
            away from the position it is pointing at, for example
            to point at the edge of a marker independent of zoom.
            Note that this shortens the arrow from the `ax` / `ay`
            vector, in contrast to `xshift` / `yshift` which moves
            everything by this amount.
        startarrowhead
            Sets the start annotation arrow head style.
        startarrowsize
            Sets the size of the start annotation arrow head,
            relative to `arrowwidth`. A value of 1 (default) gives
            a head about 3x as wide as the line.
        startstandoff
            Sets a distance, in pixels, to move the start arrowhead
            away from the position it is pointing at, for example
            to point at the edge of a marker independent of zoom.
            Note that this shortens the arrow from the `ax` / `ay`
            vector, in contrast to `xshift` / `yshift` which moves
            everything by this amount.
        templateitemname
            Used to refer to a named item in this array in the
            template. Named items from the template will be created
            even without a matching item in the input figure, but
            you can modify one by making an item with
            `templateitemname` matching its `name`, alongside your
            modifications (including `visible: false` or `enabled:
            false` to hide it). If there is no template or no
            matching item, this item will be hidden unless you
            explicitly show it with `visible: true`.
        text
            Sets the text associated with this annotation. Plotly
            uses a subset of HTML tags to do things like newline
            (`<br>`), bold (`<b></b>`), italics (`<i></i>`),
            hyperlinks (`<a href='...'></a>`). Tags `<em>`,
            `<sup>`, `<sub>`, `<s>`, `<u>`, and `<span>` are also
            supported.
        textangle
            Sets the angle at which the `text` is drawn with
            respect to the horizontal.
        valign
            Sets the vertical alignment of the `text` within the
            box. Has an effect only if an explicit height is set to
            override the text height.
        visible
            Determines whether or not this annotation is visible.
        width
            Sets an explicit width for the text box. null (default)
            lets the text set the box width. Wider text will be
            clipped. There is no automatic wrapping; use <br> to
            start a new line.
        x
            Sets the annotation's x position. If the axis `type` is
            "log", then you must take the log of your desired
            range. If the axis `type` is "date", it should be date
            strings, like date data, though Date objects and unix
            milliseconds will be accepted and converted to strings.
            If the axis `type` is "category", it should be numbers,
            using the scale where each category is assigned a
            serial number from zero in the order it appears.
        xanchor
            Sets the text box's horizontal position anchor This
            anchor binds the `x` position to the "left", "center"
            or "right" of the annotation. For example, if `x` is
            set to 1, `xref` to "paper" and `xanchor` to "right"
            then the right-most portion of the annotation lines up
            with the right-most edge of the plotting area. If
            "auto", the anchor is equivalent to "center" for data-
            referenced annotations or if there is an arrow, whereas
            for paper-referenced with no arrow, the anchor picked
            corresponds to the closest side.
        xclick
            Toggle this annotation when clicking a data point whose
            `x` value is `xclick` rather than the annotation's `x`
            value.
        xref
            Sets the annotation's x coordinate axis. If set to a x
            axis id (e.g. "x" or "x2"), the `x` position refers to
            a x coordinate. If set to "paper", the `x` position
            refers to the distance from the left of the plotting
            area in normalized coordinates where 0 (1) corresponds
            to the left (right). If set to a x axis ID followed by
            "domain" (separated by a space), the position behaves
            like for "paper", but refers to the distance in
            fractions of the domain length from the left of the
            domain of that axis: e.g., *x2 domain* refers to the
            domain of the second x  axis and a x position of 0.5
            refers to the point between the left and the right of
            the domain of the second x axis.
        xshift
            Shifts the position of the whole annotation and arrow
            to the right (positive) or left (negative) by this many
            pixels.
        y
            Sets the annotation's y position. If the axis `type` is
            "log", then you must take the log of your desired
            range. If the axis `type` is "date", it should be date
            strings, like date data, though Date objects and unix
            milliseconds will be accepted and converted to strings.
            If the axis `type` is "category", it should be numbers,
            using the scale where each category is assigned a
            serial number from zero in the order it appears.
        yanchor
            Sets the text box's vertical position anchor This
            anchor binds the `y` position to the "top", "middle" or
            "bottom" of the annotation. For example, if `y` is set
            to 1, `yref` to "paper" and `yanchor` to "top" then the
            top-most portion of the annotation lines up with the
            top-most edge of the plotting area. If "auto", the
            anchor is equivalent to "middle" for data-referenced
            annotations or if there is an arrow, whereas for paper-
            referenced with no arrow, the anchor picked corresponds
            to the closest side.
        yclick
            Toggle this annotation when clicking a data point whose
            `y` value is `yclick` rather than the annotation's `y`
            value.
        yref
            Sets the annotation's y coordinate axis. If set to a y
            axis id (e.g. "y" or "y2"), the `y` position refers to
            a y coordinate. If set to "paper", the `y` position
            refers to the distance from the bottom of the plotting
            area in normalized coordinates where 0 (1) corresponds
            to the bottom (top). If set to a y axis ID followed by
            "domain" (separated by a space), the position behaves
            like for "paper", but refers to the distance in
            fractions of the domain length from the bottom of the
            domain of that axis: e.g., *y2 domain* refers to the
            domain of the second y  axis and a y position of 0.5
            refers to the point between the bottom and the top of
            the domain of the second y axis.
        yshift
            Shifts the position of the whole annotation and arrow
            up (positive) or down (negative) by this many pixels.

        Returns
        -------
        Annotation
        """
        super().__init__("annotations")
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
The first argument to the plotly.graph_objs.layout.Annotation
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.Annotation`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("align", arg, align)
        self._set_property("arrowcolor", arg, arrowcolor)
        self._set_property("arrowhead", arg, arrowhead)
        self._set_property("arrowside", arg, arrowside)
        self._set_property("arrowsize", arg, arrowsize)
        self._set_property("arrowwidth", arg, arrowwidth)
        self._set_property("ax", arg, ax)
        self._set_property("axref", arg, axref)
        self._set_property("ay", arg, ay)
        self._set_property("ayref", arg, ayref)
        self._set_property("bgcolor", arg, bgcolor)
        self._set_property("bordercolor", arg, bordercolor)
        self._set_property("borderpad", arg, borderpad)
        self._set_property("borderwidth", arg, borderwidth)
        self._set_property("captureevents", arg, captureevents)
        self._set_property("clicktoshow", arg, clicktoshow)
        self._set_property("font", arg, font)
        self._set_property("height", arg, height)
        self._set_property("hoverlabel", arg, hoverlabel)
        self._set_property("hovertext", arg, hovertext)
        self._set_property("name", arg, name)
        self._set_property("opacity", arg, opacity)
        self._set_property("showarrow", arg, showarrow)
        self._set_property("standoff", arg, standoff)
        self._set_property("startarrowhead", arg, startarrowhead)
        self._set_property("startarrowsize", arg, startarrowsize)
        self._set_property("startstandoff", arg, startstandoff)
        self._set_property("templateitemname", arg, templateitemname)
        self._set_property("text", arg, text)
        self._set_property("textangle", arg, textangle)
        self._set_property("valign", arg, valign)
        self._set_property("visible", arg, visible)
        self._set_property("width", arg, width)
        self._set_property("x", arg, x)
        self._set_property("xanchor", arg, xanchor)
        self._set_property("xclick", arg, xclick)
        self._set_property("xref", arg, xref)
        self._set_property("xshift", arg, xshift)
        self._set_property("y", arg, y)
        self._set_property("yanchor", arg, yanchor)
        self._set_property("yclick", arg, yclick)
        self._set_property("yref", arg, yref)
        self._set_property("yshift", arg, yshift)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
