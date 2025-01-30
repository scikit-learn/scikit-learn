from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Slider(_BaseLayoutHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "layout"
    _path_str = "layout.slider"
    _valid_props = {
        "active",
        "activebgcolor",
        "bgcolor",
        "bordercolor",
        "borderwidth",
        "currentvalue",
        "font",
        "len",
        "lenmode",
        "minorticklen",
        "name",
        "pad",
        "stepdefaults",
        "steps",
        "templateitemname",
        "tickcolor",
        "ticklen",
        "tickwidth",
        "transition",
        "visible",
        "x",
        "xanchor",
        "y",
        "yanchor",
    }

    # active
    # ------
    @property
    def active(self):
        """
        Determines which button (by index starting from 0) is
        considered active.

        The 'active' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["active"]

    @active.setter
    def active(self, val):
        self["active"] = val

    # activebgcolor
    # -------------
    @property
    def activebgcolor(self):
        """
        Sets the background color of the slider grip while dragging.

        The 'activebgcolor' property is a color and may be specified as:
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
        return self["activebgcolor"]

    @activebgcolor.setter
    def activebgcolor(self, val):
        self["activebgcolor"] = val

    # bgcolor
    # -------
    @property
    def bgcolor(self):
        """
        Sets the background color of the slider.

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
        Sets the color of the border enclosing the slider.

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
        Sets the width (in px) of the border enclosing the slider.

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

    # currentvalue
    # ------------
    @property
    def currentvalue(self):
        """
        The 'currentvalue' property is an instance of Currentvalue
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.slider.Currentvalue`
          - A dict of string/value properties that will be passed
            to the Currentvalue constructor

            Supported dict properties:

                font
                    Sets the font of the current value label text.
                offset
                    The amount of space, in pixels, between the
                    current value label and the slider.
                prefix
                    When currentvalue.visible is true, this sets
                    the prefix of the label.
                suffix
                    When currentvalue.visible is true, this sets
                    the suffix of the label.
                visible
                    Shows the currently-selected value above the
                    slider.
                xanchor
                    The alignment of the value readout relative to
                    the length of the slider.

        Returns
        -------
        plotly.graph_objs.layout.slider.Currentvalue
        """
        return self["currentvalue"]

    @currentvalue.setter
    def currentvalue(self, val):
        self["currentvalue"] = val

    # font
    # ----
    @property
    def font(self):
        """
        Sets the font of the slider step labels.

        The 'font' property is an instance of Font
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.slider.Font`
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
        plotly.graph_objs.layout.slider.Font
        """
        return self["font"]

    @font.setter
    def font(self, val):
        self["font"] = val

    # len
    # ---
    @property
    def len(self):
        """
        Sets the length of the slider This measure excludes the padding
        of both ends. That is, the slider's length is this length minus
        the padding on both ends.

        The 'len' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["len"]

    @len.setter
    def len(self, val):
        self["len"] = val

    # lenmode
    # -------
    @property
    def lenmode(self):
        """
        Determines whether this slider length is set in units of plot
        "fraction" or in *pixels. Use `len` to set the value.

        The 'lenmode' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['fraction', 'pixels']

        Returns
        -------
        Any
        """
        return self["lenmode"]

    @lenmode.setter
    def lenmode(self, val):
        self["lenmode"] = val

    # minorticklen
    # ------------
    @property
    def minorticklen(self):
        """
        Sets the length in pixels of minor step tick marks

        The 'minorticklen' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["minorticklen"]

    @minorticklen.setter
    def minorticklen(self, val):
        self["minorticklen"] = val

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

    # pad
    # ---
    @property
    def pad(self):
        """
        Set the padding of the slider component along each side.

        The 'pad' property is an instance of Pad
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.slider.Pad`
          - A dict of string/value properties that will be passed
            to the Pad constructor

            Supported dict properties:

                b
                    The amount of padding (in px) along the bottom
                    of the component.
                l
                    The amount of padding (in px) on the left side
                    of the component.
                r
                    The amount of padding (in px) on the right side
                    of the component.
                t
                    The amount of padding (in px) along the top of
                    the component.

        Returns
        -------
        plotly.graph_objs.layout.slider.Pad
        """
        return self["pad"]

    @pad.setter
    def pad(self, val):
        self["pad"] = val

    # steps
    # -----
    @property
    def steps(self):
        """
        The 'steps' property is a tuple of instances of
        Step that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.slider.Step
          - A list or tuple of dicts of string/value properties that
            will be passed to the Step constructor

            Supported dict properties:

                args
                    Sets the arguments values to be passed to the
                    Plotly method set in `method` on slide.
                execute
                    When true, the API method is executed. When
                    false, all other behaviors are the same and
                    command execution is skipped. This may be
                    useful when hooking into, for example, the
                    `plotly_sliderchange` method and executing the
                    API command manually without losing the benefit
                    of the slider automatically binding to the
                    state of the plot through the specification of
                    `method` and `args`.
                label
                    Sets the text label to appear on the slider
                method
                    Sets the Plotly method to be called when the
                    slider value is changed. If the `skip` method
                    is used, the API slider will function as normal
                    but will perform no API calls and will not bind
                    automatically to state updates. This may be
                    used to create a component interface and attach
                    to slider events manually via JavaScript.
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
                value
                    Sets the value of the slider step, used to
                    refer to the step programatically. Defaults to
                    the slider label if not provided.
                visible
                    Determines whether or not this step is included
                    in the slider.

        Returns
        -------
        tuple[plotly.graph_objs.layout.slider.Step]
        """
        return self["steps"]

    @steps.setter
    def steps(self, val):
        self["steps"] = val

    # stepdefaults
    # ------------
    @property
    def stepdefaults(self):
        """
        When used in a template (as
        layout.template.layout.slider.stepdefaults), sets the default
        property values to use for elements of layout.slider.steps

        The 'stepdefaults' property is an instance of Step
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.slider.Step`
          - A dict of string/value properties that will be passed
            to the Step constructor

            Supported dict properties:

        Returns
        -------
        plotly.graph_objs.layout.slider.Step
        """
        return self["stepdefaults"]

    @stepdefaults.setter
    def stepdefaults(self, val):
        self["stepdefaults"] = val

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

    # tickcolor
    # ---------
    @property
    def tickcolor(self):
        """
        Sets the color of the border enclosing the slider.

        The 'tickcolor' property is a color and may be specified as:
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
        return self["tickcolor"]

    @tickcolor.setter
    def tickcolor(self, val):
        self["tickcolor"] = val

    # ticklen
    # -------
    @property
    def ticklen(self):
        """
        Sets the length in pixels of step tick marks

        The 'ticklen' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["ticklen"]

    @ticklen.setter
    def ticklen(self, val):
        self["ticklen"] = val

    # tickwidth
    # ---------
    @property
    def tickwidth(self):
        """
        Sets the tick width (in px).

        The 'tickwidth' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["tickwidth"]

    @tickwidth.setter
    def tickwidth(self, val):
        self["tickwidth"] = val

    # transition
    # ----------
    @property
    def transition(self):
        """
        The 'transition' property is an instance of Transition
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.slider.Transition`
          - A dict of string/value properties that will be passed
            to the Transition constructor

            Supported dict properties:

                duration
                    Sets the duration of the slider transition
                easing
                    Sets the easing function of the slider
                    transition

        Returns
        -------
        plotly.graph_objs.layout.slider.Transition
        """
        return self["transition"]

    @transition.setter
    def transition(self, val):
        self["transition"] = val

    # visible
    # -------
    @property
    def visible(self):
        """
        Determines whether or not the slider is visible.

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
        Sets the x position (in normalized coordinates) of the slider.

        The 'x' property is a number and may be specified as:
          - An int or float in the interval [-2, 3]

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
        Sets the slider's horizontal position anchor. This anchor binds
        the `x` position to the "left", "center" or "right" of the
        range selector.

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

    # y
    # -
    @property
    def y(self):
        """
        Sets the y position (in normalized coordinates) of the slider.

        The 'y' property is a number and may be specified as:
          - An int or float in the interval [-2, 3]

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
        Sets the slider's vertical position anchor This anchor binds
        the `y` position to the "top", "middle" or "bottom" of the
        range selector.

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

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        active
            Determines which button (by index starting from 0) is
            considered active.
        activebgcolor
            Sets the background color of the slider grip while
            dragging.
        bgcolor
            Sets the background color of the slider.
        bordercolor
            Sets the color of the border enclosing the slider.
        borderwidth
            Sets the width (in px) of the border enclosing the
            slider.
        currentvalue
            :class:`plotly.graph_objects.layout.slider.Currentvalue
            ` instance or dict with compatible properties
        font
            Sets the font of the slider step labels.
        len
            Sets the length of the slider This measure excludes the
            padding of both ends. That is, the slider's length is
            this length minus the padding on both ends.
        lenmode
            Determines whether this slider length is set in units
            of plot "fraction" or in *pixels. Use `len` to set the
            value.
        minorticklen
            Sets the length in pixels of minor step tick marks
        name
            When used in a template, named items are created in the
            output figure in addition to any items the figure
            already has in this array. You can modify these items
            in the output figure by making your own item with
            `templateitemname` matching this `name` alongside your
            modifications (including `visible: false` or `enabled:
            false` to hide it). Has no effect outside of a
            template.
        pad
            Set the padding of the slider component along each
            side.
        steps
            A tuple of
            :class:`plotly.graph_objects.layout.slider.Step`
            instances or dicts with compatible properties
        stepdefaults
            When used in a template (as
            layout.template.layout.slider.stepdefaults), sets the
            default property values to use for elements of
            layout.slider.steps
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
        tickcolor
            Sets the color of the border enclosing the slider.
        ticklen
            Sets the length in pixels of step tick marks
        tickwidth
            Sets the tick width (in px).
        transition
            :class:`plotly.graph_objects.layout.slider.Transition`
            instance or dict with compatible properties
        visible
            Determines whether or not the slider is visible.
        x
            Sets the x position (in normalized coordinates) of the
            slider.
        xanchor
            Sets the slider's horizontal position anchor. This
            anchor binds the `x` position to the "left", "center"
            or "right" of the range selector.
        y
            Sets the y position (in normalized coordinates) of the
            slider.
        yanchor
            Sets the slider's vertical position anchor This anchor
            binds the `y` position to the "top", "middle" or
            "bottom" of the range selector.
        """

    def __init__(
        self,
        arg=None,
        active=None,
        activebgcolor=None,
        bgcolor=None,
        bordercolor=None,
        borderwidth=None,
        currentvalue=None,
        font=None,
        len=None,
        lenmode=None,
        minorticklen=None,
        name=None,
        pad=None,
        steps=None,
        stepdefaults=None,
        templateitemname=None,
        tickcolor=None,
        ticklen=None,
        tickwidth=None,
        transition=None,
        visible=None,
        x=None,
        xanchor=None,
        y=None,
        yanchor=None,
        **kwargs,
    ):
        """
        Construct a new Slider object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of :class:`plotly.graph_objs.layout.Slider`
        active
            Determines which button (by index starting from 0) is
            considered active.
        activebgcolor
            Sets the background color of the slider grip while
            dragging.
        bgcolor
            Sets the background color of the slider.
        bordercolor
            Sets the color of the border enclosing the slider.
        borderwidth
            Sets the width (in px) of the border enclosing the
            slider.
        currentvalue
            :class:`plotly.graph_objects.layout.slider.Currentvalue
            ` instance or dict with compatible properties
        font
            Sets the font of the slider step labels.
        len
            Sets the length of the slider This measure excludes the
            padding of both ends. That is, the slider's length is
            this length minus the padding on both ends.
        lenmode
            Determines whether this slider length is set in units
            of plot "fraction" or in *pixels. Use `len` to set the
            value.
        minorticklen
            Sets the length in pixels of minor step tick marks
        name
            When used in a template, named items are created in the
            output figure in addition to any items the figure
            already has in this array. You can modify these items
            in the output figure by making your own item with
            `templateitemname` matching this `name` alongside your
            modifications (including `visible: false` or `enabled:
            false` to hide it). Has no effect outside of a
            template.
        pad
            Set the padding of the slider component along each
            side.
        steps
            A tuple of
            :class:`plotly.graph_objects.layout.slider.Step`
            instances or dicts with compatible properties
        stepdefaults
            When used in a template (as
            layout.template.layout.slider.stepdefaults), sets the
            default property values to use for elements of
            layout.slider.steps
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
        tickcolor
            Sets the color of the border enclosing the slider.
        ticklen
            Sets the length in pixels of step tick marks
        tickwidth
            Sets the tick width (in px).
        transition
            :class:`plotly.graph_objects.layout.slider.Transition`
            instance or dict with compatible properties
        visible
            Determines whether or not the slider is visible.
        x
            Sets the x position (in normalized coordinates) of the
            slider.
        xanchor
            Sets the slider's horizontal position anchor. This
            anchor binds the `x` position to the "left", "center"
            or "right" of the range selector.
        y
            Sets the y position (in normalized coordinates) of the
            slider.
        yanchor
            Sets the slider's vertical position anchor This anchor
            binds the `y` position to the "top", "middle" or
            "bottom" of the range selector.

        Returns
        -------
        Slider
        """
        super(Slider, self).__init__("sliders")

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
The first argument to the plotly.graph_objs.layout.Slider
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.Slider`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("active", None)
        _v = active if active is not None else _v
        if _v is not None:
            self["active"] = _v
        _v = arg.pop("activebgcolor", None)
        _v = activebgcolor if activebgcolor is not None else _v
        if _v is not None:
            self["activebgcolor"] = _v
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
        _v = arg.pop("currentvalue", None)
        _v = currentvalue if currentvalue is not None else _v
        if _v is not None:
            self["currentvalue"] = _v
        _v = arg.pop("font", None)
        _v = font if font is not None else _v
        if _v is not None:
            self["font"] = _v
        _v = arg.pop("len", None)
        _v = len if len is not None else _v
        if _v is not None:
            self["len"] = _v
        _v = arg.pop("lenmode", None)
        _v = lenmode if lenmode is not None else _v
        if _v is not None:
            self["lenmode"] = _v
        _v = arg.pop("minorticklen", None)
        _v = minorticklen if minorticklen is not None else _v
        if _v is not None:
            self["minorticklen"] = _v
        _v = arg.pop("name", None)
        _v = name if name is not None else _v
        if _v is not None:
            self["name"] = _v
        _v = arg.pop("pad", None)
        _v = pad if pad is not None else _v
        if _v is not None:
            self["pad"] = _v
        _v = arg.pop("steps", None)
        _v = steps if steps is not None else _v
        if _v is not None:
            self["steps"] = _v
        _v = arg.pop("stepdefaults", None)
        _v = stepdefaults if stepdefaults is not None else _v
        if _v is not None:
            self["stepdefaults"] = _v
        _v = arg.pop("templateitemname", None)
        _v = templateitemname if templateitemname is not None else _v
        if _v is not None:
            self["templateitemname"] = _v
        _v = arg.pop("tickcolor", None)
        _v = tickcolor if tickcolor is not None else _v
        if _v is not None:
            self["tickcolor"] = _v
        _v = arg.pop("ticklen", None)
        _v = ticklen if ticklen is not None else _v
        if _v is not None:
            self["ticklen"] = _v
        _v = arg.pop("tickwidth", None)
        _v = tickwidth if tickwidth is not None else _v
        if _v is not None:
            self["tickwidth"] = _v
        _v = arg.pop("transition", None)
        _v = transition if transition is not None else _v
        if _v is not None:
            self["transition"] = _v
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
        _v = arg.pop("y", None)
        _v = y if y is not None else _v
        if _v is not None:
            self["y"] = _v
        _v = arg.pop("yanchor", None)
        _v = yanchor if yanchor is not None else _v
        if _v is not None:
            self["yanchor"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
