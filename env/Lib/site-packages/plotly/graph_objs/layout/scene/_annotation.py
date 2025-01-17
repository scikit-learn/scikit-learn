from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Annotation(_BaseLayoutHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "layout.scene"
    _path_str = "layout.scene.annotation"
    _valid_props = {
        "align",
        "arrowcolor",
        "arrowhead",
        "arrowside",
        "arrowsize",
        "arrowwidth",
        "ax",
        "ay",
        "bgcolor",
        "bordercolor",
        "borderpad",
        "borderwidth",
        "captureevents",
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
        "xshift",
        "y",
        "yanchor",
        "yshift",
        "z",
    }

    # align
    # -----
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

    # arrowcolor
    # ----------
    @property
    def arrowcolor(self):
        """
        Sets the color of the annotation arrow.

        The 'arrowcolor' property is a color and may be specified as:
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
        return self["arrowcolor"]

    @arrowcolor.setter
    def arrowcolor(self, val):
        self["arrowcolor"] = val

    # arrowhead
    # ---------
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

    # arrowside
    # ---------
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

    # arrowsize
    # ---------
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

    # arrowwidth
    # ----------
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

    # ax
    # --
    @property
    def ax(self):
        """
        Sets the x component of the arrow tail about the arrow head (in
        pixels).

        The 'ax' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["ax"]

    @ax.setter
    def ax(self, val):
        self["ax"] = val

    # ay
    # --
    @property
    def ay(self):
        """
        Sets the y component of the arrow tail about the arrow head (in
        pixels).

        The 'ay' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["ay"]

    @ay.setter
    def ay(self, val):
        self["ay"] = val

    # bgcolor
    # -------
    @property
    def bgcolor(self):
        """
        Sets the background color of the annotation.

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
        Sets the color of the border enclosing the annotation `text`.

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

    # borderpad
    # ---------
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

    # borderwidth
    # -----------
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

    # captureevents
    # -------------
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

    # font
    # ----
    @property
    def font(self):
        """
        Sets the annotation text font.

        The 'font' property is an instance of Font
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.scene.annotation.Font`
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
        plotly.graph_objs.layout.scene.annotation.Font
        """
        return self["font"]

    @font.setter
    def font(self, val):
        self["font"] = val

    # height
    # ------
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

    # hoverlabel
    # ----------
    @property
    def hoverlabel(self):
        """
        The 'hoverlabel' property is an instance of Hoverlabel
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.scene.annotation.Hoverlabel`
          - A dict of string/value properties that will be passed
            to the Hoverlabel constructor

            Supported dict properties:

                bgcolor
                    Sets the background color of the hover label.
                    By default uses the annotation's `bgcolor` made
                    opaque, or white if it was transparent.
                bordercolor
                    Sets the border color of the hover label. By
                    default uses either dark grey or white, for
                    maximum contrast with `hoverlabel.bgcolor`.
                font
                    Sets the hover label text font. By default uses
                    the global hover font and size, with color from
                    `hoverlabel.bordercolor`.

        Returns
        -------
        plotly.graph_objs.layout.scene.annotation.Hoverlabel
        """
        return self["hoverlabel"]

    @hoverlabel.setter
    def hoverlabel(self, val):
        self["hoverlabel"] = val

    # hovertext
    # ---------
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

    # name
    # ----
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

    # opacity
    # -------
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

    # showarrow
    # ---------
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

    # standoff
    # --------
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

    # startarrowhead
    # --------------
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

    # startarrowsize
    # --------------
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

    # startstandoff
    # -------------
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

    # templateitemname
    # ----------------
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

    # text
    # ----
    @property
    def text(self):
        """
        Sets the text associated with this annotation. Plotly uses a
        subset of HTML tags to do things like newline (<br>), bold
        (<b></b>), italics (<i></i>), hyperlinks (<a href='...'></a>).
        Tags <em>, <sup>, <sub>, <s>, <u> <span> are also supported.

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

    # textangle
    # ---------
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

    # valign
    # ------
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

    # visible
    # -------
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

    # width
    # -----
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

    # x
    # -
    @property
    def x(self):
        """
        Sets the annotation's x position.

        The 'x' property accepts values of any type

        Returns
        -------
        Any
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

    # xshift
    # ------
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

    # y
    # -
    @property
    def y(self):
        """
        Sets the annotation's y position.

        The 'y' property accepts values of any type

        Returns
        -------
        Any
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

    # yshift
    # ------
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

    # z
    # -
    @property
    def z(self):
        """
        Sets the annotation's z position.

        The 'z' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["z"]

    @z.setter
    def z(self, val):
        self["z"] = val

    # Self properties description
    # ---------------------------
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
            head (in pixels).
        ay
            Sets the y component of the arrow tail about the arrow
            head (in pixels).
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
        font
            Sets the annotation text font.
        height
            Sets an explicit height for the text box. null
            (default) lets the text set the box height. Taller text
            will be clipped.
        hoverlabel
            :class:`plotly.graph_objects.layout.scene.annotation.Ho
            verlabel` instance or dict with compatible properties
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
            (<br>), bold (<b></b>), italics (<i></i>), hyperlinks
            (<a href='...'></a>). Tags <em>, <sup>, <sub>, <s>, <u>
            <span> are also supported.
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
            Sets the annotation's x position.
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
        xshift
            Shifts the position of the whole annotation and arrow
            to the right (positive) or left (negative) by this many
            pixels.
        y
            Sets the annotation's y position.
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
        yshift
            Shifts the position of the whole annotation and arrow
            up (positive) or down (negative) by this many pixels.
        z
            Sets the annotation's z position.
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
        ay=None,
        bgcolor=None,
        bordercolor=None,
        borderpad=None,
        borderwidth=None,
        captureevents=None,
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
        xshift=None,
        y=None,
        yanchor=None,
        yshift=None,
        z=None,
        **kwargs,
    ):
        """
        Construct a new Annotation object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.layout.scene.Annotation`
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
            head (in pixels).
        ay
            Sets the y component of the arrow tail about the arrow
            head (in pixels).
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
        font
            Sets the annotation text font.
        height
            Sets an explicit height for the text box. null
            (default) lets the text set the box height. Taller text
            will be clipped.
        hoverlabel
            :class:`plotly.graph_objects.layout.scene.annotation.Ho
            verlabel` instance or dict with compatible properties
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
            (<br>), bold (<b></b>), italics (<i></i>), hyperlinks
            (<a href='...'></a>). Tags <em>, <sup>, <sub>, <s>, <u>
            <span> are also supported.
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
            Sets the annotation's x position.
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
        xshift
            Shifts the position of the whole annotation and arrow
            to the right (positive) or left (negative) by this many
            pixels.
        y
            Sets the annotation's y position.
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
        yshift
            Shifts the position of the whole annotation and arrow
            up (positive) or down (negative) by this many pixels.
        z
            Sets the annotation's z position.

        Returns
        -------
        Annotation
        """
        super(Annotation, self).__init__("annotations")

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
The first argument to the plotly.graph_objs.layout.scene.Annotation
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.scene.Annotation`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("align", None)
        _v = align if align is not None else _v
        if _v is not None:
            self["align"] = _v
        _v = arg.pop("arrowcolor", None)
        _v = arrowcolor if arrowcolor is not None else _v
        if _v is not None:
            self["arrowcolor"] = _v
        _v = arg.pop("arrowhead", None)
        _v = arrowhead if arrowhead is not None else _v
        if _v is not None:
            self["arrowhead"] = _v
        _v = arg.pop("arrowside", None)
        _v = arrowside if arrowside is not None else _v
        if _v is not None:
            self["arrowside"] = _v
        _v = arg.pop("arrowsize", None)
        _v = arrowsize if arrowsize is not None else _v
        if _v is not None:
            self["arrowsize"] = _v
        _v = arg.pop("arrowwidth", None)
        _v = arrowwidth if arrowwidth is not None else _v
        if _v is not None:
            self["arrowwidth"] = _v
        _v = arg.pop("ax", None)
        _v = ax if ax is not None else _v
        if _v is not None:
            self["ax"] = _v
        _v = arg.pop("ay", None)
        _v = ay if ay is not None else _v
        if _v is not None:
            self["ay"] = _v
        _v = arg.pop("bgcolor", None)
        _v = bgcolor if bgcolor is not None else _v
        if _v is not None:
            self["bgcolor"] = _v
        _v = arg.pop("bordercolor", None)
        _v = bordercolor if bordercolor is not None else _v
        if _v is not None:
            self["bordercolor"] = _v
        _v = arg.pop("borderpad", None)
        _v = borderpad if borderpad is not None else _v
        if _v is not None:
            self["borderpad"] = _v
        _v = arg.pop("borderwidth", None)
        _v = borderwidth if borderwidth is not None else _v
        if _v is not None:
            self["borderwidth"] = _v
        _v = arg.pop("captureevents", None)
        _v = captureevents if captureevents is not None else _v
        if _v is not None:
            self["captureevents"] = _v
        _v = arg.pop("font", None)
        _v = font if font is not None else _v
        if _v is not None:
            self["font"] = _v
        _v = arg.pop("height", None)
        _v = height if height is not None else _v
        if _v is not None:
            self["height"] = _v
        _v = arg.pop("hoverlabel", None)
        _v = hoverlabel if hoverlabel is not None else _v
        if _v is not None:
            self["hoverlabel"] = _v
        _v = arg.pop("hovertext", None)
        _v = hovertext if hovertext is not None else _v
        if _v is not None:
            self["hovertext"] = _v
        _v = arg.pop("name", None)
        _v = name if name is not None else _v
        if _v is not None:
            self["name"] = _v
        _v = arg.pop("opacity", None)
        _v = opacity if opacity is not None else _v
        if _v is not None:
            self["opacity"] = _v
        _v = arg.pop("showarrow", None)
        _v = showarrow if showarrow is not None else _v
        if _v is not None:
            self["showarrow"] = _v
        _v = arg.pop("standoff", None)
        _v = standoff if standoff is not None else _v
        if _v is not None:
            self["standoff"] = _v
        _v = arg.pop("startarrowhead", None)
        _v = startarrowhead if startarrowhead is not None else _v
        if _v is not None:
            self["startarrowhead"] = _v
        _v = arg.pop("startarrowsize", None)
        _v = startarrowsize if startarrowsize is not None else _v
        if _v is not None:
            self["startarrowsize"] = _v
        _v = arg.pop("startstandoff", None)
        _v = startstandoff if startstandoff is not None else _v
        if _v is not None:
            self["startstandoff"] = _v
        _v = arg.pop("templateitemname", None)
        _v = templateitemname if templateitemname is not None else _v
        if _v is not None:
            self["templateitemname"] = _v
        _v = arg.pop("text", None)
        _v = text if text is not None else _v
        if _v is not None:
            self["text"] = _v
        _v = arg.pop("textangle", None)
        _v = textangle if textangle is not None else _v
        if _v is not None:
            self["textangle"] = _v
        _v = arg.pop("valign", None)
        _v = valign if valign is not None else _v
        if _v is not None:
            self["valign"] = _v
        _v = arg.pop("visible", None)
        _v = visible if visible is not None else _v
        if _v is not None:
            self["visible"] = _v
        _v = arg.pop("width", None)
        _v = width if width is not None else _v
        if _v is not None:
            self["width"] = _v
        _v = arg.pop("x", None)
        _v = x if x is not None else _v
        if _v is not None:
            self["x"] = _v
        _v = arg.pop("xanchor", None)
        _v = xanchor if xanchor is not None else _v
        if _v is not None:
            self["xanchor"] = _v
        _v = arg.pop("xshift", None)
        _v = xshift if xshift is not None else _v
        if _v is not None:
            self["xshift"] = _v
        _v = arg.pop("y", None)
        _v = y if y is not None else _v
        if _v is not None:
            self["y"] = _v
        _v = arg.pop("yanchor", None)
        _v = yanchor if yanchor is not None else _v
        if _v is not None:
            self["yanchor"] = _v
        _v = arg.pop("yshift", None)
        _v = yshift if yshift is not None else _v
        if _v is not None:
            self["yshift"] = _v
        _v = arg.pop("z", None)
        _v = z if z is not None else _v
        if _v is not None:
            self["z"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
