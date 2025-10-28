#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Hoverlabel(_BaseTraceHierarchyType):
    _parent_path_str = "densitymap"
    _path_str = "densitymap.hoverlabel"
    _valid_props = {
        "align",
        "alignsrc",
        "bgcolor",
        "bgcolorsrc",
        "bordercolor",
        "bordercolorsrc",
        "font",
        "namelength",
        "namelengthsrc",
        "showarrow",
    }

    @property
    def align(self):
        """
        Sets the horizontal alignment of the text content within hover
        label box. Has an effect only if the hover label text spans
        more two or more lines

        The 'align' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['left', 'right', 'auto']
          - A tuple, list, or one-dimensional numpy array of the above

        Returns
        -------
        Any|numpy.ndarray
        """
        return self["align"]

    @align.setter
    def align(self, val):
        self["align"] = val

    @property
    def alignsrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `align`.

        The 'alignsrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["alignsrc"]

    @alignsrc.setter
    def alignsrc(self, val):
        self["alignsrc"] = val

    @property
    def bgcolor(self):
        """
        Sets the background color of the hover labels for this trace

        The 'bgcolor' property is a color and may be specified as:
          - A hex string (e.g. '#ff0000')
          - An rgb/rgba string (e.g. 'rgb(255,0,0)')
          - An hsl/hsla string (e.g. 'hsl(0,100%,50%)')
          - An hsv/hsva string (e.g. 'hsv(0,100%,100%)')
          - A named CSS color: see https://plotly.com/python/css-colors/ for a list
          - A list or array of any of the above

        Returns
        -------
        str|numpy.ndarray
        """
        return self["bgcolor"]

    @bgcolor.setter
    def bgcolor(self, val):
        self["bgcolor"] = val

    @property
    def bgcolorsrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `bgcolor`.

        The 'bgcolorsrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["bgcolorsrc"]

    @bgcolorsrc.setter
    def bgcolorsrc(self, val):
        self["bgcolorsrc"] = val

    @property
    def bordercolor(self):
        """
        Sets the border color of the hover labels for this trace.

        The 'bordercolor' property is a color and may be specified as:
          - A hex string (e.g. '#ff0000')
          - An rgb/rgba string (e.g. 'rgb(255,0,0)')
          - An hsl/hsla string (e.g. 'hsl(0,100%,50%)')
          - An hsv/hsva string (e.g. 'hsv(0,100%,100%)')
          - A named CSS color: see https://plotly.com/python/css-colors/ for a list
          - A list or array of any of the above

        Returns
        -------
        str|numpy.ndarray
        """
        return self["bordercolor"]

    @bordercolor.setter
    def bordercolor(self, val):
        self["bordercolor"] = val

    @property
    def bordercolorsrc(self):
        """
        Sets the source reference on Chart Studio Cloud for
        `bordercolor`.

        The 'bordercolorsrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["bordercolorsrc"]

    @bordercolorsrc.setter
    def bordercolorsrc(self, val):
        self["bordercolorsrc"] = val

    @property
    def font(self):
        """
        Sets the font used in hover labels.

        The 'font' property is an instance of Font
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.densitymap.hoverlabel.Font`
          - A dict of string/value properties that will be passed
            to the Font constructor

        Returns
        -------
        plotly.graph_objs.densitymap.hoverlabel.Font
        """
        return self["font"]

    @font.setter
    def font(self, val):
        self["font"] = val

    @property
    def namelength(self):
        """
        Sets the default length (in number of characters) of the trace
        name in the hover labels for all traces. -1 shows the whole
        name regardless of length. 0-3 shows the first 0-3 characters,
        and an integer >3 will show the whole name if it is less than
        that many characters, but if it is longer, will truncate to
        `namelength - 3` characters and add an ellipsis.

        The 'namelength' property is a integer and may be specified as:
          - An int (or float that will be cast to an int)
            in the interval [-1, 9223372036854775807]
          - A tuple, list, or one-dimensional numpy array of the above

        Returns
        -------
        int|numpy.ndarray
        """
        return self["namelength"]

    @namelength.setter
    def namelength(self, val):
        self["namelength"] = val

    @property
    def namelengthsrc(self):
        """
        Sets the source reference on Chart Studio Cloud for
        `namelength`.

        The 'namelengthsrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["namelengthsrc"]

    @namelengthsrc.setter
    def namelengthsrc(self, val):
        self["namelengthsrc"] = val

    @property
    def showarrow(self):
        """
        Sets whether or not to show the hover label arrow/triangle
        pointing to the data point.

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
    def _prop_descriptions(self):
        return """\
        align
            Sets the horizontal alignment of the text content
            within hover label box. Has an effect only if the hover
            label text spans more two or more lines
        alignsrc
            Sets the source reference on Chart Studio Cloud for
            `align`.
        bgcolor
            Sets the background color of the hover labels for this
            trace
        bgcolorsrc
            Sets the source reference on Chart Studio Cloud for
            `bgcolor`.
        bordercolor
            Sets the border color of the hover labels for this
            trace.
        bordercolorsrc
            Sets the source reference on Chart Studio Cloud for
            `bordercolor`.
        font
            Sets the font used in hover labels.
        namelength
            Sets the default length (in number of characters) of
            the trace name in the hover labels for all traces. -1
            shows the whole name regardless of length. 0-3 shows
            the first 0-3 characters, and an integer >3 will show
            the whole name if it is less than that many characters,
            but if it is longer, will truncate to `namelength - 3`
            characters and add an ellipsis.
        namelengthsrc
            Sets the source reference on Chart Studio Cloud for
            `namelength`.
        showarrow
            Sets whether or not to show the hover label
            arrow/triangle pointing to the data point.
        """

    def __init__(
        self,
        arg=None,
        align=None,
        alignsrc=None,
        bgcolor=None,
        bgcolorsrc=None,
        bordercolor=None,
        bordercolorsrc=None,
        font=None,
        namelength=None,
        namelengthsrc=None,
        showarrow=None,
        **kwargs,
    ):
        """
        Construct a new Hoverlabel object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.densitymap.Hoverlabel`
        align
            Sets the horizontal alignment of the text content
            within hover label box. Has an effect only if the hover
            label text spans more two or more lines
        alignsrc
            Sets the source reference on Chart Studio Cloud for
            `align`.
        bgcolor
            Sets the background color of the hover labels for this
            trace
        bgcolorsrc
            Sets the source reference on Chart Studio Cloud for
            `bgcolor`.
        bordercolor
            Sets the border color of the hover labels for this
            trace.
        bordercolorsrc
            Sets the source reference on Chart Studio Cloud for
            `bordercolor`.
        font
            Sets the font used in hover labels.
        namelength
            Sets the default length (in number of characters) of
            the trace name in the hover labels for all traces. -1
            shows the whole name regardless of length. 0-3 shows
            the first 0-3 characters, and an integer >3 will show
            the whole name if it is less than that many characters,
            but if it is longer, will truncate to `namelength - 3`
            characters and add an ellipsis.
        namelengthsrc
            Sets the source reference on Chart Studio Cloud for
            `namelength`.
        showarrow
            Sets whether or not to show the hover label
            arrow/triangle pointing to the data point.

        Returns
        -------
        Hoverlabel
        """
        super().__init__("hoverlabel")
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
The first argument to the plotly.graph_objs.densitymap.Hoverlabel
constructor must be a dict or
an instance of :class:`plotly.graph_objs.densitymap.Hoverlabel`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("align", arg, align)
        self._set_property("alignsrc", arg, alignsrc)
        self._set_property("bgcolor", arg, bgcolor)
        self._set_property("bgcolorsrc", arg, bgcolorsrc)
        self._set_property("bordercolor", arg, bordercolor)
        self._set_property("bordercolorsrc", arg, bordercolorsrc)
        self._set_property("font", arg, font)
        self._set_property("namelength", arg, namelength)
        self._set_property("namelengthsrc", arg, namelengthsrc)
        self._set_property("showarrow", arg, showarrow)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
