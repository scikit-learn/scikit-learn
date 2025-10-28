#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Updatemenu(_BaseLayoutHierarchyType):
    _parent_path_str = "layout"
    _path_str = "layout.updatemenu"
    _valid_props = {
        "active",
        "bgcolor",
        "bordercolor",
        "borderwidth",
        "buttondefaults",
        "buttons",
        "direction",
        "font",
        "name",
        "pad",
        "showactive",
        "templateitemname",
        "type",
        "visible",
        "x",
        "xanchor",
        "y",
        "yanchor",
    }

    @property
    def active(self):
        """
        Determines which button (by index starting from 0) is
        considered active.

        The 'active' property is a integer and may be specified as:
          - An int (or float that will be cast to an int)
            in the interval [-1, 9223372036854775807]

        Returns
        -------
        int
        """
        return self["active"]

    @active.setter
    def active(self, val):
        self["active"] = val

    @property
    def bgcolor(self):
        """
        Sets the background color of the update menu buttons.

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
        Sets the color of the border enclosing the update menu.

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
        Sets the width (in px) of the border enclosing the update menu.

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
        The 'buttons' property is a tuple of instances of
        Button that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.updatemenu.Button
          - A list or tuple of dicts of string/value properties that
            will be passed to the Button constructor

        Returns
        -------
        tuple[plotly.graph_objs.layout.updatemenu.Button]
        """
        return self["buttons"]

    @buttons.setter
    def buttons(self, val):
        self["buttons"] = val

    @property
    def buttondefaults(self):
        """
        When used in a template (as
        layout.template.layout.updatemenu.buttondefaults), sets the
        default property values to use for elements of
        layout.updatemenu.buttons

        The 'buttondefaults' property is an instance of Button
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.updatemenu.Button`
          - A dict of string/value properties that will be passed
            to the Button constructor

        Returns
        -------
        plotly.graph_objs.layout.updatemenu.Button
        """
        return self["buttondefaults"]

    @buttondefaults.setter
    def buttondefaults(self, val):
        self["buttondefaults"] = val

    @property
    def direction(self):
        """
        Determines the direction in which the buttons are laid out,
        whether in a dropdown menu or a row/column of buttons. For
        `left` and `up`, the buttons will still appear in left-to-right
        or top-to-bottom order respectively.

        The 'direction' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['left', 'right', 'up', 'down']

        Returns
        -------
        Any
        """
        return self["direction"]

    @direction.setter
    def direction(self, val):
        self["direction"] = val

    @property
    def font(self):
        """
        Sets the font of the update menu button text.

        The 'font' property is an instance of Font
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.updatemenu.Font`
          - A dict of string/value properties that will be passed
            to the Font constructor

        Returns
        -------
        plotly.graph_objs.layout.updatemenu.Font
        """
        return self["font"]

    @font.setter
    def font(self, val):
        self["font"] = val

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
    def pad(self):
        """
        Sets the padding around the buttons or dropdown menu.

        The 'pad' property is an instance of Pad
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.updatemenu.Pad`
          - A dict of string/value properties that will be passed
            to the Pad constructor

        Returns
        -------
        plotly.graph_objs.layout.updatemenu.Pad
        """
        return self["pad"]

    @pad.setter
    def pad(self, val):
        self["pad"] = val

    @property
    def showactive(self):
        """
        Highlights active dropdown item or active button if true.

        The 'showactive' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["showactive"]

    @showactive.setter
    def showactive(self, val):
        self["showactive"] = val

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
    def type(self):
        """
        Determines whether the buttons are accessible via a dropdown
        menu or whether the buttons are stacked horizontally or
        vertically

        The 'type' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['dropdown', 'buttons']

        Returns
        -------
        Any
        """
        return self["type"]

    @type.setter
    def type(self, val):
        self["type"] = val

    @property
    def visible(self):
        """
        Determines whether or not the update menu is visible.

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
        Sets the x position (in normalized coordinates) of the update
        menu.

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
        Sets the update menu's horizontal position anchor. This anchor
        binds the `x` position to the "left", "center" or "right" of
        the range selector.

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
        Sets the y position (in normalized coordinates) of the update
        menu.

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
        Sets the update menu's vertical position anchor This anchor
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
        active
            Determines which button (by index starting from 0) is
            considered active.
        bgcolor
            Sets the background color of the update menu buttons.
        bordercolor
            Sets the color of the border enclosing the update menu.
        borderwidth
            Sets the width (in px) of the border enclosing the
            update menu.
        buttons
            A tuple of
            :class:`plotly.graph_objects.layout.updatemenu.Button`
            instances or dicts with compatible properties
        buttondefaults
            When used in a template (as
            layout.template.layout.updatemenu.buttondefaults), sets
            the default property values to use for elements of
            layout.updatemenu.buttons
        direction
            Determines the direction in which the buttons are laid
            out, whether in a dropdown menu or a row/column of
            buttons. For `left` and `up`, the buttons will still
            appear in left-to-right or top-to-bottom order
            respectively.
        font
            Sets the font of the update menu button text.
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
            Sets the padding around the buttons or dropdown menu.
        showactive
            Highlights active dropdown item or active button if
            true.
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
        type
            Determines whether the buttons are accessible via a
            dropdown menu or whether the buttons are stacked
            horizontally or vertically
        visible
            Determines whether or not the update menu is visible.
        x
            Sets the x position (in normalized coordinates) of the
            update menu.
        xanchor
            Sets the update menu's horizontal position anchor. This
            anchor binds the `x` position to the "left", "center"
            or "right" of the range selector.
        y
            Sets the y position (in normalized coordinates) of the
            update menu.
        yanchor
            Sets the update menu's vertical position anchor This
            anchor binds the `y` position to the "top", "middle" or
            "bottom" of the range selector.
        """

    def __init__(
        self,
        arg=None,
        active=None,
        bgcolor=None,
        bordercolor=None,
        borderwidth=None,
        buttons=None,
        buttondefaults=None,
        direction=None,
        font=None,
        name=None,
        pad=None,
        showactive=None,
        templateitemname=None,
        type=None,
        visible=None,
        x=None,
        xanchor=None,
        y=None,
        yanchor=None,
        **kwargs,
    ):
        """
        Construct a new Updatemenu object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.layout.Updatemenu`
        active
            Determines which button (by index starting from 0) is
            considered active.
        bgcolor
            Sets the background color of the update menu buttons.
        bordercolor
            Sets the color of the border enclosing the update menu.
        borderwidth
            Sets the width (in px) of the border enclosing the
            update menu.
        buttons
            A tuple of
            :class:`plotly.graph_objects.layout.updatemenu.Button`
            instances or dicts with compatible properties
        buttondefaults
            When used in a template (as
            layout.template.layout.updatemenu.buttondefaults), sets
            the default property values to use for elements of
            layout.updatemenu.buttons
        direction
            Determines the direction in which the buttons are laid
            out, whether in a dropdown menu or a row/column of
            buttons. For `left` and `up`, the buttons will still
            appear in left-to-right or top-to-bottom order
            respectively.
        font
            Sets the font of the update menu button text.
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
            Sets the padding around the buttons or dropdown menu.
        showactive
            Highlights active dropdown item or active button if
            true.
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
        type
            Determines whether the buttons are accessible via a
            dropdown menu or whether the buttons are stacked
            horizontally or vertically
        visible
            Determines whether or not the update menu is visible.
        x
            Sets the x position (in normalized coordinates) of the
            update menu.
        xanchor
            Sets the update menu's horizontal position anchor. This
            anchor binds the `x` position to the "left", "center"
            or "right" of the range selector.
        y
            Sets the y position (in normalized coordinates) of the
            update menu.
        yanchor
            Sets the update menu's vertical position anchor This
            anchor binds the `y` position to the "top", "middle" or
            "bottom" of the range selector.

        Returns
        -------
        Updatemenu
        """
        super().__init__("updatemenus")
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
The first argument to the plotly.graph_objs.layout.Updatemenu
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.Updatemenu`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("active", arg, active)
        self._set_property("bgcolor", arg, bgcolor)
        self._set_property("bordercolor", arg, bordercolor)
        self._set_property("borderwidth", arg, borderwidth)
        self._set_property("buttons", arg, buttons)
        self._set_property("buttondefaults", arg, buttondefaults)
        self._set_property("direction", arg, direction)
        self._set_property("font", arg, font)
        self._set_property("name", arg, name)
        self._set_property("pad", arg, pad)
        self._set_property("showactive", arg, showactive)
        self._set_property("templateitemname", arg, templateitemname)
        self._set_property("type", arg, type)
        self._set_property("visible", arg, visible)
        self._set_property("x", arg, x)
        self._set_property("xanchor", arg, xanchor)
        self._set_property("y", arg, y)
        self._set_property("yanchor", arg, yanchor)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
