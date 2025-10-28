#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Rangeselector(_BaseLayoutHierarchyType):
    _parent_path_str = "layout.xaxis"
    _path_str = "layout.xaxis.rangeselector"
    _valid_props = {
        "activecolor",
        "bgcolor",
        "bordercolor",
        "borderwidth",
        "buttondefaults",
        "buttons",
        "font",
        "visible",
        "x",
        "xanchor",
        "y",
        "yanchor",
    }

    @property
    def activecolor(self):
        """
        Sets the background color of the active range selector button.

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
    def bgcolor(self):
        """
        Sets the background color of the range selector buttons.

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
        Sets the color of the border enclosing the range selector.

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
        Sets the width (in px) of the border enclosing the range
        selector.

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
    def buttons(self):
        """
        Sets the specifications for each buttons. By default, a range
        selector comes with no buttons.

        The 'buttons' property is a tuple of instances of
        Button that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.xaxis.rangeselector.Button
          - A list or tuple of dicts of string/value properties that
            will be passed to the Button constructor

        Returns
        -------
        tuple[plotly.graph_objs.layout.xaxis.rangeselector.Button]
        """
        return self["buttons"]

    @buttons.setter
    def buttons(self, val):
        self["buttons"] = val

    @property
    def buttondefaults(self):
        """
        When used in a template (as
        layout.template.layout.xaxis.rangeselector.buttondefaults),
        sets the default property values to use for elements of
        layout.xaxis.rangeselector.buttons

        The 'buttondefaults' property is an instance of Button
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.xaxis.rangeselector.Button`
          - A dict of string/value properties that will be passed
            to the Button constructor

        Returns
        -------
        plotly.graph_objs.layout.xaxis.rangeselector.Button
        """
        return self["buttondefaults"]

    @buttondefaults.setter
    def buttondefaults(self, val):
        self["buttondefaults"] = val

    @property
    def font(self):
        """
        Sets the font of the range selector button text.

        The 'font' property is an instance of Font
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.xaxis.rangeselector.Font`
          - A dict of string/value properties that will be passed
            to the Font constructor

        Returns
        -------
        plotly.graph_objs.layout.xaxis.rangeselector.Font
        """
        return self["font"]

    @font.setter
    def font(self, val):
        self["font"] = val

    @property
    def visible(self):
        """
        Determines whether or not this range selector is visible. Note
        that range selectors are only available for x axes of `type`
        set to or auto-typed to "date".

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
        Sets the x position (in normalized coordinates) of the range
        selector.

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

    @property
    def xanchor(self):
        """
        Sets the range selector's horizontal position anchor. This
        anchor binds the `x` position to the "left", "center" or
        "right" of the range selector.

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
    def y(self):
        """
        Sets the y position (in normalized coordinates) of the range
        selector.

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

    @property
    def yanchor(self):
        """
        Sets the range selector's vertical position anchor This anchor
        binds the `y` position to the "top", "middle" or "bottom" of
        the range selector.

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
    def _prop_descriptions(self):
        return """\
        activecolor
            Sets the background color of the active range selector
            button.
        bgcolor
            Sets the background color of the range selector
            buttons.
        bordercolor
            Sets the color of the border enclosing the range
            selector.
        borderwidth
            Sets the width (in px) of the border enclosing the
            range selector.
        buttons
            Sets the specifications for each buttons. By default, a
            range selector comes with no buttons.
        buttondefaults
            When used in a template (as layout.template.layout.xaxi
            s.rangeselector.buttondefaults), sets the default
            property values to use for elements of
            layout.xaxis.rangeselector.buttons
        font
            Sets the font of the range selector button text.
        visible
            Determines whether or not this range selector is
            visible. Note that range selectors are only available
            for x axes of `type` set to or auto-typed to "date".
        x
            Sets the x position (in normalized coordinates) of the
            range selector.
        xanchor
            Sets the range selector's horizontal position anchor.
            This anchor binds the `x` position to the "left",
            "center" or "right" of the range selector.
        y
            Sets the y position (in normalized coordinates) of the
            range selector.
        yanchor
            Sets the range selector's vertical position anchor This
            anchor binds the `y` position to the "top", "middle" or
            "bottom" of the range selector.
        """

    def __init__(
        self,
        arg=None,
        activecolor=None,
        bgcolor=None,
        bordercolor=None,
        borderwidth=None,
        buttons=None,
        buttondefaults=None,
        font=None,
        visible=None,
        x=None,
        xanchor=None,
        y=None,
        yanchor=None,
        **kwargs,
    ):
        """
        Construct a new Rangeselector object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.layout.xaxis.Rangeselector`
        activecolor
            Sets the background color of the active range selector
            button.
        bgcolor
            Sets the background color of the range selector
            buttons.
        bordercolor
            Sets the color of the border enclosing the range
            selector.
        borderwidth
            Sets the width (in px) of the border enclosing the
            range selector.
        buttons
            Sets the specifications for each buttons. By default, a
            range selector comes with no buttons.
        buttondefaults
            When used in a template (as layout.template.layout.xaxi
            s.rangeselector.buttondefaults), sets the default
            property values to use for elements of
            layout.xaxis.rangeselector.buttons
        font
            Sets the font of the range selector button text.
        visible
            Determines whether or not this range selector is
            visible. Note that range selectors are only available
            for x axes of `type` set to or auto-typed to "date".
        x
            Sets the x position (in normalized coordinates) of the
            range selector.
        xanchor
            Sets the range selector's horizontal position anchor.
            This anchor binds the `x` position to the "left",
            "center" or "right" of the range selector.
        y
            Sets the y position (in normalized coordinates) of the
            range selector.
        yanchor
            Sets the range selector's vertical position anchor This
            anchor binds the `y` position to the "top", "middle" or
            "bottom" of the range selector.

        Returns
        -------
        Rangeselector
        """
        super().__init__("rangeselector")
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
The first argument to the plotly.graph_objs.layout.xaxis.Rangeselector
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.xaxis.Rangeselector`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("activecolor", arg, activecolor)
        self._set_property("bgcolor", arg, bgcolor)
        self._set_property("bordercolor", arg, bordercolor)
        self._set_property("borderwidth", arg, borderwidth)
        self._set_property("buttons", arg, buttons)
        self._set_property("buttondefaults", arg, buttondefaults)
        self._set_property("font", arg, font)
        self._set_property("visible", arg, visible)
        self._set_property("x", arg, x)
        self._set_property("xanchor", arg, xanchor)
        self._set_property("y", arg, y)
        self._set_property("yanchor", arg, yanchor)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
