from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Title(_BaseLayoutHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "layout"
    _path_str = "layout.title"
    _valid_props = {
        "automargin",
        "font",
        "pad",
        "subtitle",
        "text",
        "x",
        "xanchor",
        "xref",
        "y",
        "yanchor",
        "yref",
    }

    # automargin
    # ----------
    @property
    def automargin(self):
        """
        Determines whether the title can automatically push the figure
        margins. If `yref='paper'` then the margin will expand to
        ensure that the title doesn’t overlap with the edges of the
        container. If `yref='container'` then the margins will ensure
        that the title doesn’t overlap with the plot area, tick labels,
        and axis titles. If `automargin=true` and the margins need to
        be expanded, then y will be set to a default 1 and yanchor will
        be set to an appropriate default to ensure that minimal margin
        space is needed. Note that when `yref='paper'`, only 1 or 0 are
        allowed y values. Invalid values will be reset to the default
        1.

        The 'automargin' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["automargin"]

    @automargin.setter
    def automargin(self, val):
        self["automargin"] = val

    # font
    # ----
    @property
    def font(self):
        """
        Sets the title font. Note that the title's font used to be
        customized by the now deprecated `titlefont` attribute.

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
        plotly.graph_objs.layout.title.Font
        """
        return self["font"]

    @font.setter
    def font(self, val):
        self["font"] = val

    # pad
    # ---
    @property
    def pad(self):
        """
        Sets the padding of the title. Each padding value only applies
        when the corresponding `xanchor`/`yanchor` value is set
        accordingly. E.g. for left padding to take effect, `xanchor`
        must be set to "left". The same rule applies if
        `xanchor`/`yanchor` is determined automatically. Padding is
        muted if the respective anchor value is "middle*/*center".

        The 'pad' property is an instance of Pad
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.title.Pad`
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
        plotly.graph_objs.layout.title.Pad
        """
        return self["pad"]

    @pad.setter
    def pad(self, val):
        self["pad"] = val

    # subtitle
    # --------
    @property
    def subtitle(self):
        """
        The 'subtitle' property is an instance of Subtitle
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.title.Subtitle`
          - A dict of string/value properties that will be passed
            to the Subtitle constructor

            Supported dict properties:

                font
                    Sets the subtitle font.
                text
                    Sets the plot's subtitle.

        Returns
        -------
        plotly.graph_objs.layout.title.Subtitle
        """
        return self["subtitle"]

    @subtitle.setter
    def subtitle(self, val):
        self["subtitle"] = val

    # text
    # ----
    @property
    def text(self):
        """
        Sets the plot's title. Note that before the existence of
        `title.text`, the title's contents used to be defined as the
        `title` attribute itself. This behavior has been deprecated.

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

    # x
    # -
    @property
    def x(self):
        """
        Sets the x position with respect to `xref` in normalized
        coordinates from 0 (left) to 1 (right).

        The 'x' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

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
        Sets the title's horizontal alignment with respect to its x
        position. "left" means that the title starts at x, "right"
        means that the title ends at x and "center" means that the
        title's center is at x. "auto" divides `xref` by three and
        calculates the `xanchor` value automatically based on the value
        of `x`.

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
        Sets the y position with respect to `yref` in normalized
        coordinates from 0 (bottom) to 1 (top). "auto" places the
        baseline of the title onto the vertical center of the top
        margin.

        The 'y' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

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
        Sets the title's vertical alignment with respect to its y
        position. "top" means that the title's cap line is at y,
        "bottom" means that the title's baseline is at y and "middle"
        means that the title's midline is at y. "auto" divides `yref`
        by three and calculates the `yanchor` value automatically based
        on the value of `y`.

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
        automargin
            Determines whether the title can automatically push the
            figure margins. If `yref='paper'` then the margin will
            expand to ensure that the title doesn’t overlap with
            the edges of the container. If `yref='container'` then
            the margins will ensure that the title doesn’t overlap
            with the plot area, tick labels, and axis titles. If
            `automargin=true` and the margins need to be expanded,
            then y will be set to a default 1 and yanchor will be
            set to an appropriate default to ensure that minimal
            margin space is needed. Note that when `yref='paper'`,
            only 1 or 0 are allowed y values. Invalid values will
            be reset to the default 1.
        font
            Sets the title font. Note that the title's font used to
            be customized by the now deprecated `titlefont`
            attribute.
        pad
            Sets the padding of the title. Each padding value only
            applies when the corresponding `xanchor`/`yanchor`
            value is set accordingly. E.g. for left padding to take
            effect, `xanchor` must be set to "left". The same rule
            applies if `xanchor`/`yanchor` is determined
            automatically. Padding is muted if the respective
            anchor value is "middle*/*center".
        subtitle
            :class:`plotly.graph_objects.layout.title.Subtitle`
            instance or dict with compatible properties
        text
            Sets the plot's title. Note that before the existence
            of `title.text`, the title's contents used to be
            defined as the `title` attribute itself. This behavior
            has been deprecated.
        x
            Sets the x position with respect to `xref` in
            normalized coordinates from 0 (left) to 1 (right).
        xanchor
            Sets the title's horizontal alignment with respect to
            its x position. "left" means that the title starts at
            x, "right" means that the title ends at x and "center"
            means that the title's center is at x. "auto" divides
            `xref` by three and calculates the `xanchor` value
            automatically based on the value of `x`.
        xref
            Sets the container `x` refers to. "container" spans the
            entire `width` of the plot. "paper" refers to the width
            of the plotting area only.
        y
            Sets the y position with respect to `yref` in
            normalized coordinates from 0 (bottom) to 1 (top).
            "auto" places the baseline of the title onto the
            vertical center of the top margin.
        yanchor
            Sets the title's vertical alignment with respect to its
            y position. "top" means that the title's cap line is at
            y, "bottom" means that the title's baseline is at y and
            "middle" means that the title's midline is at y. "auto"
            divides `yref` by three and calculates the `yanchor`
            value automatically based on the value of `y`.
        yref
            Sets the container `y` refers to. "container" spans the
            entire `height` of the plot. "paper" refers to the
            height of the plotting area only.
        """

    def __init__(
        self,
        arg=None,
        automargin=None,
        font=None,
        pad=None,
        subtitle=None,
        text=None,
        x=None,
        xanchor=None,
        xref=None,
        y=None,
        yanchor=None,
        yref=None,
        **kwargs,
    ):
        """
        Construct a new Title object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of :class:`plotly.graph_objs.layout.Title`
        automargin
            Determines whether the title can automatically push the
            figure margins. If `yref='paper'` then the margin will
            expand to ensure that the title doesn’t overlap with
            the edges of the container. If `yref='container'` then
            the margins will ensure that the title doesn’t overlap
            with the plot area, tick labels, and axis titles. If
            `automargin=true` and the margins need to be expanded,
            then y will be set to a default 1 and yanchor will be
            set to an appropriate default to ensure that minimal
            margin space is needed. Note that when `yref='paper'`,
            only 1 or 0 are allowed y values. Invalid values will
            be reset to the default 1.
        font
            Sets the title font. Note that the title's font used to
            be customized by the now deprecated `titlefont`
            attribute.
        pad
            Sets the padding of the title. Each padding value only
            applies when the corresponding `xanchor`/`yanchor`
            value is set accordingly. E.g. for left padding to take
            effect, `xanchor` must be set to "left". The same rule
            applies if `xanchor`/`yanchor` is determined
            automatically. Padding is muted if the respective
            anchor value is "middle*/*center".
        subtitle
            :class:`plotly.graph_objects.layout.title.Subtitle`
            instance or dict with compatible properties
        text
            Sets the plot's title. Note that before the existence
            of `title.text`, the title's contents used to be
            defined as the `title` attribute itself. This behavior
            has been deprecated.
        x
            Sets the x position with respect to `xref` in
            normalized coordinates from 0 (left) to 1 (right).
        xanchor
            Sets the title's horizontal alignment with respect to
            its x position. "left" means that the title starts at
            x, "right" means that the title ends at x and "center"
            means that the title's center is at x. "auto" divides
            `xref` by three and calculates the `xanchor` value
            automatically based on the value of `x`.
        xref
            Sets the container `x` refers to. "container" spans the
            entire `width` of the plot. "paper" refers to the width
            of the plotting area only.
        y
            Sets the y position with respect to `yref` in
            normalized coordinates from 0 (bottom) to 1 (top).
            "auto" places the baseline of the title onto the
            vertical center of the top margin.
        yanchor
            Sets the title's vertical alignment with respect to its
            y position. "top" means that the title's cap line is at
            y, "bottom" means that the title's baseline is at y and
            "middle" means that the title's midline is at y. "auto"
            divides `yref` by three and calculates the `yanchor`
            value automatically based on the value of `y`.
        yref
            Sets the container `y` refers to. "container" spans the
            entire `height` of the plot. "paper" refers to the
            height of the plotting area only.

        Returns
        -------
        Title
        """
        super(Title, self).__init__("title")

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
The first argument to the plotly.graph_objs.layout.Title
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.Title`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("automargin", None)
        _v = automargin if automargin is not None else _v
        if _v is not None:
            self["automargin"] = _v
        _v = arg.pop("font", None)
        _v = font if font is not None else _v
        if _v is not None:
            self["font"] = _v
        _v = arg.pop("pad", None)
        _v = pad if pad is not None else _v
        if _v is not None:
            self["pad"] = _v
        _v = arg.pop("subtitle", None)
        _v = subtitle if subtitle is not None else _v
        if _v is not None:
            self["subtitle"] = _v
        _v = arg.pop("text", None)
        _v = text if text is not None else _v
        if _v is not None:
            self["text"] = _v
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
