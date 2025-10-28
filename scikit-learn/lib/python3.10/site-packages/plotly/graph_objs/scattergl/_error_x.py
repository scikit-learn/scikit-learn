#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class ErrorX(_BaseTraceHierarchyType):
    _parent_path_str = "scattergl"
    _path_str = "scattergl.error_x"
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

    @property
    def color(self):
        """
        Sets the stroke color of the error bars.

        The 'color' property is a color and may be specified as:
          - A hex string (e.g. '#ff0000')
          - An rgb/rgba string (e.g. 'rgb(255,0,0)')
          - An hsl/hsla string (e.g. 'hsl(0,100%,50%)')
          - An hsv/hsva string (e.g. 'hsv(0,100%,100%)')
          - A named CSS color: see https://plotly.com/python/css-colors/ for a list

        Returns
        -------
        str
        """
        return self["color"]

    @color.setter
    def color(self, val):
        self["color"] = val

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

    @property
    def type(self):
        """
        Determines the rule used to generate the error bars. If
        "constant", the bar lengths are of a constant value. Set this
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
            Sets the stroke color of the error bars.
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
            "constant", the bar lengths are of a constant value.
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
            :class:`plotly.graph_objs.scattergl.ErrorX`
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
            Sets the stroke color of the error bars.
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
            "constant", the bar lengths are of a constant value.
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
        super().__init__("error_x")
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
The first argument to the plotly.graph_objs.scattergl.ErrorX
constructor must be a dict or
an instance of :class:`plotly.graph_objs.scattergl.ErrorX`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("array", arg, array)
        self._set_property("arrayminus", arg, arrayminus)
        self._set_property("arrayminussrc", arg, arrayminussrc)
        self._set_property("arraysrc", arg, arraysrc)
        self._set_property("color", arg, color)
        self._set_property("copy_ystyle", arg, copy_ystyle)
        self._set_property("symmetric", arg, symmetric)
        self._set_property("thickness", arg, thickness)
        self._set_property("traceref", arg, traceref)
        self._set_property("tracerefminus", arg, tracerefminus)
        self._set_property("type", arg, type)
        self._set_property("value", arg, value)
        self._set_property("valueminus", arg, valueminus)
        self._set_property("visible", arg, visible)
        self._set_property("width", arg, width)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
