#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Title(_BaseLayoutHierarchyType):
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

    @property
    def font(self):
        """
        Sets the title font.

        The 'font' property is an instance of Font
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.title.Font`
          - A dict of string/value properties that will be passed
            to the Font constructor

        Returns
        -------
        plotly.graph_objs.layout.title.Font
        """
        return self["font"]

    @font.setter
    def font(self, val):
        self["font"] = val

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

        Returns
        -------
        plotly.graph_objs.layout.title.Pad
        """
        return self["pad"]

    @pad.setter
    def pad(self, val):
        self["pad"] = val

    @property
    def subtitle(self):
        """
        The 'subtitle' property is an instance of Subtitle
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.title.Subtitle`
          - A dict of string/value properties that will be passed
            to the Subtitle constructor

        Returns
        -------
        plotly.graph_objs.layout.title.Subtitle
        """
        return self["subtitle"]

    @subtitle.setter
    def subtitle(self, val):
        self["subtitle"] = val

    @property
    def text(self):
        """
        Sets the plot's title.

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
            Sets the title font.
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
            Sets the plot's title.
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
            Sets the title font.
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
            Sets the plot's title.
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
        super().__init__("title")
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
The first argument to the plotly.graph_objs.layout.Title
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.Title`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("automargin", arg, automargin)
        self._set_property("font", arg, font)
        self._set_property("pad", arg, pad)
        self._set_property("subtitle", arg, subtitle)
        self._set_property("text", arg, text)
        self._set_property("x", arg, x)
        self._set_property("xanchor", arg, xanchor)
        self._set_property("xref", arg, xref)
        self._set_property("y", arg, y)
        self._set_property("yanchor", arg, yanchor)
        self._set_property("yref", arg, yref)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
