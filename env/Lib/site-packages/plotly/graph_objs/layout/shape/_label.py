from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Label(_BaseLayoutHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "layout.shape"
    _path_str = "layout.shape.label"
    _valid_props = {
        "font",
        "padding",
        "text",
        "textangle",
        "textposition",
        "texttemplate",
        "xanchor",
        "yanchor",
    }

    # font
    # ----
    @property
    def font(self):
        """
        Sets the shape label text font.

        The 'font' property is an instance of Font
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.shape.label.Font`
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
        plotly.graph_objs.layout.shape.label.Font
        """
        return self["font"]

    @font.setter
    def font(self, val):
        self["font"] = val

    # padding
    # -------
    @property
    def padding(self):
        """
        Sets padding (in px) between edge of label and edge of shape.

        The 'padding' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["padding"]

    @padding.setter
    def padding(self, val):
        self["padding"] = val

    # text
    # ----
    @property
    def text(self):
        """
        Sets the text to display with shape. It is also used for legend
        item if `name` is not provided.

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
        Sets the angle at which the label text is drawn with respect to
        the horizontal. For lines, angle "auto" is the same angle as
        the line. For all other shapes, angle "auto" is horizontal.

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

    # textposition
    # ------------
    @property
    def textposition(self):
        """
        Sets the position of the label text relative to the shape.
        Supported values for rectangles, circles and paths are *top
        left*, *top center*, *top right*, *middle left*, *middle
        center*, *middle right*, *bottom left*, *bottom center*, and
        *bottom right*. Supported values for lines are "start",
        "middle", and "end". Default: *middle center* for rectangles,
        circles, and paths; "middle" for lines.

        The 'textposition' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['top left', 'top center', 'top right', 'middle left',
                'middle center', 'middle right', 'bottom left', 'bottom
                center', 'bottom right', 'start', 'middle', 'end']

        Returns
        -------
        Any
        """
        return self["textposition"]

    @textposition.setter
    def textposition(self, val):
        self["textposition"] = val

    # texttemplate
    # ------------
    @property
    def texttemplate(self):
        """
        Template string used for rendering the shape's label. Note that
        this will override `text`. Variables are inserted using
        %{variable}, for example "x0: %{x0}". Numbers are formatted
        using d3-format's syntax %{variable:d3-format}, for example
        "Price: %{x0:$.2f}". See
        https://github.com/d3/d3-format/tree/v1.4.5#d3-format for
        details on the formatting syntax. Dates are formatted using
        d3-time-format's syntax %{variable|d3-time-format}, for example
        "Day: %{x0|%m %b %Y}". See https://github.com/d3/d3-time-
        format/tree/v2.2.3#locale_format for details on the date
        formatting syntax. A single multiplication or division
        operation may be applied to numeric variables, and combined
        with d3 number formatting, for example "Length in cm:
        %{x0*2.54}", "%{slope*60:.1f} meters per second." For log axes,
        variable values are given in log units. For date axes, x/y
        coordinate variables and center variables use datetimes, while
        all other variable values use values in ms. Finally, the
        template string has access to variables `x0`, `x1`, `y0`, `y1`,
        `slope`, `dx`, `dy`, `width`, `height`, `length`, `xcenter` and
        `ycenter`.

        The 'texttemplate' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["texttemplate"]

    @texttemplate.setter
    def texttemplate(self, val):
        self["texttemplate"] = val

    # xanchor
    # -------
    @property
    def xanchor(self):
        """
        Sets the label's horizontal position anchor This anchor binds
        the specified `textposition` to the "left", "center" or "right"
        of the label text. For example, if `textposition` is set to
        *top right* and `xanchor` to "right" then the right-most
        portion of the label text lines up with the right-most edge of
        the shape.

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

    # yanchor
    # -------
    @property
    def yanchor(self):
        """
        Sets the label's vertical position anchor This anchor binds the
        specified `textposition` to the "top", "middle" or "bottom" of
        the label text. For example, if `textposition` is set to *top
        right* and `yanchor` to "top" then the top-most portion of the
        label text lines up with the top-most edge of the shape.

        The 'yanchor' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['top', 'middle', 'bottom']

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
        font
            Sets the shape label text font.
        padding
            Sets padding (in px) between edge of label and edge of
            shape.
        text
            Sets the text to display with shape. It is also used
            for legend item if `name` is not provided.
        textangle
            Sets the angle at which the label text is drawn with
            respect to the horizontal. For lines, angle "auto" is
            the same angle as the line. For all other shapes, angle
            "auto" is horizontal.
        textposition
            Sets the position of the label text relative to the
            shape. Supported values for rectangles, circles and
            paths are *top left*, *top center*, *top right*,
            *middle left*, *middle center*, *middle right*, *bottom
            left*, *bottom center*, and *bottom right*. Supported
            values for lines are "start", "middle", and "end".
            Default: *middle center* for rectangles, circles, and
            paths; "middle" for lines.
        texttemplate
            Template string used for rendering the shape's label.
            Note that this will override `text`. Variables are
            inserted using %{variable}, for example "x0: %{x0}".
            Numbers are formatted using d3-format's syntax
            %{variable:d3-format}, for example "Price: %{x0:$.2f}".
            See
            https://github.com/d3/d3-format/tree/v1.4.5#d3-format
            for details on the formatting syntax. Dates are
            formatted using d3-time-format's syntax
            %{variable|d3-time-format}, for example "Day: %{x0|%m
            %b %Y}". See https://github.com/d3/d3-time-
            format/tree/v2.2.3#locale_format for details on the
            date formatting syntax. A single multiplication or
            division operation may be applied to numeric variables,
            and combined with d3 number formatting, for example
            "Length in cm: %{x0*2.54}", "%{slope*60:.1f} meters per
            second." For log axes, variable values are given in log
            units. For date axes, x/y coordinate variables and
            center variables use datetimes, while all other
            variable values use values in ms. Finally, the template
            string has access to variables `x0`, `x1`, `y0`, `y1`,
            `slope`, `dx`, `dy`, `width`, `height`, `length`,
            `xcenter` and `ycenter`.
        xanchor
            Sets the label's horizontal position anchor This anchor
            binds the specified `textposition` to the "left",
            "center" or "right" of the label text. For example, if
            `textposition` is set to *top right* and `xanchor` to
            "right" then the right-most portion of the label text
            lines up with the right-most edge of the shape.
        yanchor
            Sets the label's vertical position anchor This anchor
            binds the specified `textposition` to the "top",
            "middle" or "bottom" of the label text. For example, if
            `textposition` is set to *top right* and `yanchor` to
            "top" then the top-most portion of the label text lines
            up with the top-most edge of the shape.
        """

    def __init__(
        self,
        arg=None,
        font=None,
        padding=None,
        text=None,
        textangle=None,
        textposition=None,
        texttemplate=None,
        xanchor=None,
        yanchor=None,
        **kwargs,
    ):
        """
        Construct a new Label object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.layout.shape.Label`
        font
            Sets the shape label text font.
        padding
            Sets padding (in px) between edge of label and edge of
            shape.
        text
            Sets the text to display with shape. It is also used
            for legend item if `name` is not provided.
        textangle
            Sets the angle at which the label text is drawn with
            respect to the horizontal. For lines, angle "auto" is
            the same angle as the line. For all other shapes, angle
            "auto" is horizontal.
        textposition
            Sets the position of the label text relative to the
            shape. Supported values for rectangles, circles and
            paths are *top left*, *top center*, *top right*,
            *middle left*, *middle center*, *middle right*, *bottom
            left*, *bottom center*, and *bottom right*. Supported
            values for lines are "start", "middle", and "end".
            Default: *middle center* for rectangles, circles, and
            paths; "middle" for lines.
        texttemplate
            Template string used for rendering the shape's label.
            Note that this will override `text`. Variables are
            inserted using %{variable}, for example "x0: %{x0}".
            Numbers are formatted using d3-format's syntax
            %{variable:d3-format}, for example "Price: %{x0:$.2f}".
            See
            https://github.com/d3/d3-format/tree/v1.4.5#d3-format
            for details on the formatting syntax. Dates are
            formatted using d3-time-format's syntax
            %{variable|d3-time-format}, for example "Day: %{x0|%m
            %b %Y}". See https://github.com/d3/d3-time-
            format/tree/v2.2.3#locale_format for details on the
            date formatting syntax. A single multiplication or
            division operation may be applied to numeric variables,
            and combined with d3 number formatting, for example
            "Length in cm: %{x0*2.54}", "%{slope*60:.1f} meters per
            second." For log axes, variable values are given in log
            units. For date axes, x/y coordinate variables and
            center variables use datetimes, while all other
            variable values use values in ms. Finally, the template
            string has access to variables `x0`, `x1`, `y0`, `y1`,
            `slope`, `dx`, `dy`, `width`, `height`, `length`,
            `xcenter` and `ycenter`.
        xanchor
            Sets the label's horizontal position anchor This anchor
            binds the specified `textposition` to the "left",
            "center" or "right" of the label text. For example, if
            `textposition` is set to *top right* and `xanchor` to
            "right" then the right-most portion of the label text
            lines up with the right-most edge of the shape.
        yanchor
            Sets the label's vertical position anchor This anchor
            binds the specified `textposition` to the "top",
            "middle" or "bottom" of the label text. For example, if
            `textposition` is set to *top right* and `yanchor` to
            "top" then the top-most portion of the label text lines
            up with the top-most edge of the shape.

        Returns
        -------
        Label
        """
        super(Label, self).__init__("label")

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
The first argument to the plotly.graph_objs.layout.shape.Label
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.shape.Label`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("font", None)
        _v = font if font is not None else _v
        if _v is not None:
            self["font"] = _v
        _v = arg.pop("padding", None)
        _v = padding if padding is not None else _v
        if _v is not None:
            self["padding"] = _v
        _v = arg.pop("text", None)
        _v = text if text is not None else _v
        if _v is not None:
            self["text"] = _v
        _v = arg.pop("textangle", None)
        _v = textangle if textangle is not None else _v
        if _v is not None:
            self["textangle"] = _v
        _v = arg.pop("textposition", None)
        _v = textposition if textposition is not None else _v
        if _v is not None:
            self["textposition"] = _v
        _v = arg.pop("texttemplate", None)
        _v = texttemplate if texttemplate is not None else _v
        if _v is not None:
            self["texttemplate"] = _v
        _v = arg.pop("xanchor", None)
        _v = xanchor if xanchor is not None else _v
        if _v is not None:
            self["xanchor"] = _v
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
