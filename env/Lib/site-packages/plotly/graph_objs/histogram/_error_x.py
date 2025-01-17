from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class ErrorX(_BaseTraceHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "histogram"
    _path_str = "histogram.error_x"
    _valid_props = {
        "array",
        "arrayminus",
        "arrayminussrc",
        "arraysrc",
        "color",
        "copy_ystyle",
        "symmetric",
        "thickness",
        "traceref",
        "tracerefminus",
        "type",
        "value",
        "valueminus",
        "visible",
        "width",
    }

    # array
    # -----
    @property
    def array(self):
        """
        Sets the data corresponding the length of each error bar.
        Values are plotted relative to the underlying data.

        The 'array' property is an array that may be specified as a tuple,
        list, numpy array, or pandas Series

        Returns
        -------
        numpy.ndarray
        """
        return self["array"]

    @array.setter
    def array(self, val):
        self["array"] = val

    # arrayminus
    # ----------
    @property
    def arrayminus(self):
        """
        Sets the data corresponding the length of each error bar in the
        bottom (left) direction for vertical (horizontal) bars Values
        are plotted relative to the underlying data.

        The 'arrayminus' property is an array that may be specified as a tuple,
        list, numpy array, or pandas Series

        Returns
        -------
        numpy.ndarray
        """
        return self["arrayminus"]

    @arrayminus.setter
    def arrayminus(self, val):
        self["arrayminus"] = val

    # arrayminussrc
    # -------------
    @property
    def arrayminussrc(self):
        """
        Sets the source reference on Chart Studio Cloud for
        `arrayminus`.

        The 'arrayminussrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["arrayminussrc"]

    @arrayminussrc.setter
    def arrayminussrc(self, val):
        self["arrayminussrc"] = val

    # arraysrc
    # --------
    @property
    def arraysrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `array`.

        The 'arraysrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["arraysrc"]

    @arraysrc.setter
    def arraysrc(self, val):
        self["arraysrc"] = val

    # color
    # -----
    @property
    def color(self):
        """
        Sets the stoke color of the error bars.

        The 'color' property is a color and may be specified as:
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
        return self["color"]

    @color.setter
    def color(self, val):
        self["color"] = val

    # copy_ystyle
    # -----------
    @property
    def copy_ystyle(self):
        """
        The 'copy_ystyle' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["copy_ystyle"]

    @copy_ystyle.setter
    def copy_ystyle(self, val):
        self["copy_ystyle"] = val

    # symmetric
    # ---------
    @property
    def symmetric(self):
        """
        Determines whether or not the error bars have the same length
        in both direction (top/bottom for vertical bars, left/right for
        horizontal bars.

        The 'symmetric' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["symmetric"]

    @symmetric.setter
    def symmetric(self, val):
        self["symmetric"] = val

    # thickness
    # ---------
    @property
    def thickness(self):
        """
        Sets the thickness (in px) of the error bars.

        The 'thickness' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["thickness"]

    @thickness.setter
    def thickness(self, val):
        self["thickness"] = val

    # traceref
    # --------
    @property
    def traceref(self):
        """
        The 'traceref' property is a integer and may be specified as:
          - An int (or float that will be cast to an int)
            in the interval [0, 9223372036854775807]

        Returns
        -------
        int
        """
        return self["traceref"]

    @traceref.setter
    def traceref(self, val):
        self["traceref"] = val

    # tracerefminus
    # -------------
    @property
    def tracerefminus(self):
        """
        The 'tracerefminus' property is a integer and may be specified as:
          - An int (or float that will be cast to an int)
            in the interval [0, 9223372036854775807]

        Returns
        -------
        int
        """
        return self["tracerefminus"]

    @tracerefminus.setter
    def tracerefminus(self, val):
        self["tracerefminus"] = val

    # type
    # ----
    @property
    def type(self):
        """
        Determines the rule used to generate the error bars. If
        *constant`, the bar lengths are of a constant value. Set this
        constant in `value`. If "percent", the bar lengths correspond
        to a percentage of underlying data. Set this percentage in
        `value`. If "sqrt", the bar lengths correspond to the square of
        the underlying data. If "data", the bar lengths are set with
        data set `array`.

        The 'type' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['percent', 'constant', 'sqrt', 'data']

        Returns
        -------
        Any
        """
        return self["type"]

    @type.setter
    def type(self, val):
        self["type"] = val

    # value
    # -----
    @property
    def value(self):
        """
        Sets the value of either the percentage (if `type` is set to
        "percent") or the constant (if `type` is set to "constant")
        corresponding to the lengths of the error bars.

        The 'value' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["value"]

    @value.setter
    def value(self, val):
        self["value"] = val

    # valueminus
    # ----------
    @property
    def valueminus(self):
        """
        Sets the value of either the percentage (if `type` is set to
        "percent") or the constant (if `type` is set to "constant")
        corresponding to the lengths of the error bars in the bottom
        (left) direction for vertical (horizontal) bars

        The 'valueminus' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["valueminus"]

    @valueminus.setter
    def valueminus(self, val):
        self["valueminus"] = val

    # visible
    # -------
    @property
    def visible(self):
        """
        Determines whether or not this set of error bars is visible.

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
        Sets the width (in px) of the cross-bar at both ends of the
        error bars.

        The 'width' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["width"]

    @width.setter
    def width(self, val):
        self["width"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        array
            Sets the data corresponding the length of each error
            bar. Values are plotted relative to the underlying
            data.
        arrayminus
            Sets the data corresponding the length of each error
            bar in the bottom (left) direction for vertical
            (horizontal) bars Values are plotted relative to the
            underlying data.
        arrayminussrc
            Sets the source reference on Chart Studio Cloud for
            `arrayminus`.
        arraysrc
            Sets the source reference on Chart Studio Cloud for
            `array`.
        color
            Sets the stoke color of the error bars.
        copy_ystyle

        symmetric
            Determines whether or not the error bars have the same
            length in both direction (top/bottom for vertical bars,
            left/right for horizontal bars.
        thickness
            Sets the thickness (in px) of the error bars.
        traceref

        tracerefminus

        type
            Determines the rule used to generate the error bars. If
            *constant`, the bar lengths are of a constant value.
            Set this constant in `value`. If "percent", the bar
            lengths correspond to a percentage of underlying data.
            Set this percentage in `value`. If "sqrt", the bar
            lengths correspond to the square of the underlying
            data. If "data", the bar lengths are set with data set
            `array`.
        value
            Sets the value of either the percentage (if `type` is
            set to "percent") or the constant (if `type` is set to
            "constant") corresponding to the lengths of the error
            bars.
        valueminus
            Sets the value of either the percentage (if `type` is
            set to "percent") or the constant (if `type` is set to
            "constant") corresponding to the lengths of the error
            bars in the bottom (left) direction for vertical
            (horizontal) bars
        visible
            Determines whether or not this set of error bars is
            visible.
        width
            Sets the width (in px) of the cross-bar at both ends of
            the error bars.
        """

    def __init__(
        self,
        arg=None,
        array=None,
        arrayminus=None,
        arrayminussrc=None,
        arraysrc=None,
        color=None,
        copy_ystyle=None,
        symmetric=None,
        thickness=None,
        traceref=None,
        tracerefminus=None,
        type=None,
        value=None,
        valueminus=None,
        visible=None,
        width=None,
        **kwargs,
    ):
        """
        Construct a new ErrorX object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.histogram.ErrorX`
        array
            Sets the data corresponding the length of each error
            bar. Values are plotted relative to the underlying
            data.
        arrayminus
            Sets the data corresponding the length of each error
            bar in the bottom (left) direction for vertical
            (horizontal) bars Values are plotted relative to the
            underlying data.
        arrayminussrc
            Sets the source reference on Chart Studio Cloud for
            `arrayminus`.
        arraysrc
            Sets the source reference on Chart Studio Cloud for
            `array`.
        color
            Sets the stoke color of the error bars.
        copy_ystyle

        symmetric
            Determines whether or not the error bars have the same
            length in both direction (top/bottom for vertical bars,
            left/right for horizontal bars.
        thickness
            Sets the thickness (in px) of the error bars.
        traceref

        tracerefminus

        type
            Determines the rule used to generate the error bars. If
            *constant`, the bar lengths are of a constant value.
            Set this constant in `value`. If "percent", the bar
            lengths correspond to a percentage of underlying data.
            Set this percentage in `value`. If "sqrt", the bar
            lengths correspond to the square of the underlying
            data. If "data", the bar lengths are set with data set
            `array`.
        value
            Sets the value of either the percentage (if `type` is
            set to "percent") or the constant (if `type` is set to
            "constant") corresponding to the lengths of the error
            bars.
        valueminus
            Sets the value of either the percentage (if `type` is
            set to "percent") or the constant (if `type` is set to
            "constant") corresponding to the lengths of the error
            bars in the bottom (left) direction for vertical
            (horizontal) bars
        visible
            Determines whether or not this set of error bars is
            visible.
        width
            Sets the width (in px) of the cross-bar at both ends of
            the error bars.

        Returns
        -------
        ErrorX
        """
        super(ErrorX, self).__init__("error_x")

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
The first argument to the plotly.graph_objs.histogram.ErrorX
constructor must be a dict or
an instance of :class:`plotly.graph_objs.histogram.ErrorX`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("array", None)
        _v = array if array is not None else _v
        if _v is not None:
            self["array"] = _v
        _v = arg.pop("arrayminus", None)
        _v = arrayminus if arrayminus is not None else _v
        if _v is not None:
            self["arrayminus"] = _v
        _v = arg.pop("arrayminussrc", None)
        _v = arrayminussrc if arrayminussrc is not None else _v
        if _v is not None:
            self["arrayminussrc"] = _v
        _v = arg.pop("arraysrc", None)
        _v = arraysrc if arraysrc is not None else _v
        if _v is not None:
            self["arraysrc"] = _v
        _v = arg.pop("color", None)
        _v = color if color is not None else _v
        if _v is not None:
            self["color"] = _v
        _v = arg.pop("copy_ystyle", None)
        _v = copy_ystyle if copy_ystyle is not None else _v
        if _v is not None:
            self["copy_ystyle"] = _v
        _v = arg.pop("symmetric", None)
        _v = symmetric if symmetric is not None else _v
        if _v is not None:
            self["symmetric"] = _v
        _v = arg.pop("thickness", None)
        _v = thickness if thickness is not None else _v
        if _v is not None:
            self["thickness"] = _v
        _v = arg.pop("traceref", None)
        _v = traceref if traceref is not None else _v
        if _v is not None:
            self["traceref"] = _v
        _v = arg.pop("tracerefminus", None)
        _v = tracerefminus if tracerefminus is not None else _v
        if _v is not None:
            self["tracerefminus"] = _v
        _v = arg.pop("type", None)
        _v = type if type is not None else _v
        if _v is not None:
            self["type"] = _v
        _v = arg.pop("value", None)
        _v = value if value is not None else _v
        if _v is not None:
            self["value"] = _v
        _v = arg.pop("valueminus", None)
        _v = valueminus if valueminus is not None else _v
        if _v is not None:
            self["valueminus"] = _v
        _v = arg.pop("visible", None)
        _v = visible if visible is not None else _v
        if _v is not None:
            self["visible"] = _v
        _v = arg.pop("width", None)
        _v = width if width is not None else _v
        if _v is not None:
            self["width"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
