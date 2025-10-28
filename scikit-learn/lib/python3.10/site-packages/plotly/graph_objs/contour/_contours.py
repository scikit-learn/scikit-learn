#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Contours(_BaseTraceHierarchyType):
    _parent_path_str = "contour"
    _path_str = "contour.contours"
    _valid_props = {
        "coloring",
        "end",
        "labelfont",
        "labelformat",
        "operation",
        "showlabels",
        "showlines",
        "size",
        "start",
        "type",
        "value",
    }

    @property
    def coloring(self):
        """
        Determines the coloring method showing the contour values. If
        "fill", coloring is done evenly between each contour level If
        "heatmap", a heatmap gradient coloring is applied between each
        contour level. If "lines", coloring is done on the contour
        lines. If "none", no coloring is applied on this trace.

        The 'coloring' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['fill', 'heatmap', 'lines', 'none']

        Returns
        -------
        Any
        """
        return self["coloring"]

    @coloring.setter
    def coloring(self, val):
        self["coloring"] = val

    @property
    def end(self):
        """
        Sets the end contour level value. Must be more than
        `contours.start`

        The 'end' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["end"]

    @end.setter
    def end(self, val):
        self["end"] = val

    @property
    def labelfont(self):
        """
        Sets the font used for labeling the contour levels. The default
        color comes from the lines, if shown. The default family and
        size come from `layout.font`.

        The 'labelfont' property is an instance of Labelfont
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.contour.contours.Labelfont`
          - A dict of string/value properties that will be passed
            to the Labelfont constructor

        Returns
        -------
        plotly.graph_objs.contour.contours.Labelfont
        """
        return self["labelfont"]

    @labelfont.setter
    def labelfont(self, val):
        self["labelfont"] = val

    @property
    def labelformat(self):
        """
        Sets the contour label formatting rule using d3 formatting
        mini-languages which are very similar to those in Python. For
        numbers, see:
        https://github.com/d3/d3-format/tree/v1.4.5#d3-format.

        The 'labelformat' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["labelformat"]

    @labelformat.setter
    def labelformat(self, val):
        self["labelformat"] = val

    @property
    def operation(self):
        """
        Sets the constraint operation. "=" keeps regions equal to
        `value` "<" and "<=" keep regions less than `value` ">" and
        ">=" keep regions greater than `value` "[]", "()", "[)", and
        "(]" keep regions inside `value[0]` to `value[1]` "][", ")(",
        "](", ")[" keep regions outside `value[0]` to value[1]` Open
        vs. closed intervals make no difference to constraint display,
        but all versions are allowed for consistency with filter
        transforms.

        The 'operation' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['=', '<', '>=', '>', '<=', '[]', '()', '[)', '(]', '][',
                ')(', '](', ')[']

        Returns
        -------
        Any
        """
        return self["operation"]

    @operation.setter
    def operation(self, val):
        self["operation"] = val

    @property
    def showlabels(self):
        """
        Determines whether to label the contour lines with their
        values.

        The 'showlabels' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["showlabels"]

    @showlabels.setter
    def showlabels(self, val):
        self["showlabels"] = val

    @property
    def showlines(self):
        """
        Determines whether or not the contour lines are drawn. Has an
        effect only if `contours.coloring` is set to "fill".

        The 'showlines' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["showlines"]

    @showlines.setter
    def showlines(self, val):
        self["showlines"] = val

    @property
    def size(self):
        """
        Sets the step between each contour level. Must be positive.

        The 'size' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["size"]

    @size.setter
    def size(self, val):
        self["size"] = val

    @property
    def start(self):
        """
        Sets the starting contour level value. Must be less than
        `contours.end`

        The 'start' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["start"]

    @start.setter
    def start(self, val):
        self["start"] = val

    @property
    def type(self):
        """
        If `levels`, the data is represented as a contour plot with
        multiple levels displayed. If `constraint`, the data is
        represented as constraints with the invalid region shaded as
        specified by the `operation` and `value` parameters.

        The 'type' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['levels', 'constraint']

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
        Sets the value or values of the constraint boundary. When
        `operation` is set to one of the comparison values
        (`=,<,>=,>,<=`) "value" is expected to be a number. When
        `operation` is set to one of the interval values
        (`[],(),[),(],][,)(,](,)[`) "value" is expected to be an array
        of two numbers where the first is the lower bound and the
        second is the upper bound.

        The 'value' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["value"]

    @value.setter
    def value(self, val):
        self["value"] = val

    @property
    def _prop_descriptions(self):
        return """\
        coloring
            Determines the coloring method showing the contour
            values. If "fill", coloring is done evenly between each
            contour level If "heatmap", a heatmap gradient coloring
            is applied between each contour level. If "lines",
            coloring is done on the contour lines. If "none", no
            coloring is applied on this trace.
        end
            Sets the end contour level value. Must be more than
            `contours.start`
        labelfont
            Sets the font used for labeling the contour levels. The
            default color comes from the lines, if shown. The
            default family and size come from `layout.font`.
        labelformat
            Sets the contour label formatting rule using d3
            formatting mini-languages which are very similar to
            those in Python. For numbers, see:
            https://github.com/d3/d3-format/tree/v1.4.5#d3-format.
        operation
            Sets the constraint operation. "=" keeps regions equal
            to `value` "<" and "<=" keep regions less than `value`
            ">" and ">=" keep regions greater than `value` "[]",
            "()", "[)", and "(]" keep regions inside `value[0]` to
            `value[1]` "][", ")(", "](", ")[" keep regions outside
            `value[0]` to value[1]` Open vs. closed intervals make
            no difference to constraint display, but all versions
            are allowed for consistency with filter transforms.
        showlabels
            Determines whether to label the contour lines with
            their values.
        showlines
            Determines whether or not the contour lines are drawn.
            Has an effect only if `contours.coloring` is set to
            "fill".
        size
            Sets the step between each contour level. Must be
            positive.
        start
            Sets the starting contour level value. Must be less
            than `contours.end`
        type
            If `levels`, the data is represented as a contour plot
            with multiple levels displayed. If `constraint`, the
            data is represented as constraints with the invalid
            region shaded as specified by the `operation` and
            `value` parameters.
        value
            Sets the value or values of the constraint boundary.
            When `operation` is set to one of the comparison values
            (`=,<,>=,>,<=`) "value" is expected to be a number.
            When `operation` is set to one of the interval values
            (`[],(),[),(],][,)(,](,)[`) "value" is expected to be
            an array of two numbers where the first is the lower
            bound and the second is the upper bound.
        """

    def __init__(
        self,
        arg=None,
        coloring=None,
        end=None,
        labelfont=None,
        labelformat=None,
        operation=None,
        showlabels=None,
        showlines=None,
        size=None,
        start=None,
        type=None,
        value=None,
        **kwargs,
    ):
        """
        Construct a new Contours object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.contour.Contours`
        coloring
            Determines the coloring method showing the contour
            values. If "fill", coloring is done evenly between each
            contour level If "heatmap", a heatmap gradient coloring
            is applied between each contour level. If "lines",
            coloring is done on the contour lines. If "none", no
            coloring is applied on this trace.
        end
            Sets the end contour level value. Must be more than
            `contours.start`
        labelfont
            Sets the font used for labeling the contour levels. The
            default color comes from the lines, if shown. The
            default family and size come from `layout.font`.
        labelformat
            Sets the contour label formatting rule using d3
            formatting mini-languages which are very similar to
            those in Python. For numbers, see:
            https://github.com/d3/d3-format/tree/v1.4.5#d3-format.
        operation
            Sets the constraint operation. "=" keeps regions equal
            to `value` "<" and "<=" keep regions less than `value`
            ">" and ">=" keep regions greater than `value` "[]",
            "()", "[)", and "(]" keep regions inside `value[0]` to
            `value[1]` "][", ")(", "](", ")[" keep regions outside
            `value[0]` to value[1]` Open vs. closed intervals make
            no difference to constraint display, but all versions
            are allowed for consistency with filter transforms.
        showlabels
            Determines whether to label the contour lines with
            their values.
        showlines
            Determines whether or not the contour lines are drawn.
            Has an effect only if `contours.coloring` is set to
            "fill".
        size
            Sets the step between each contour level. Must be
            positive.
        start
            Sets the starting contour level value. Must be less
            than `contours.end`
        type
            If `levels`, the data is represented as a contour plot
            with multiple levels displayed. If `constraint`, the
            data is represented as constraints with the invalid
            region shaded as specified by the `operation` and
            `value` parameters.
        value
            Sets the value or values of the constraint boundary.
            When `operation` is set to one of the comparison values
            (`=,<,>=,>,<=`) "value" is expected to be a number.
            When `operation` is set to one of the interval values
            (`[],(),[),(],][,)(,](,)[`) "value" is expected to be
            an array of two numbers where the first is the lower
            bound and the second is the upper bound.

        Returns
        -------
        Contours
        """
        super().__init__("contours")
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
The first argument to the plotly.graph_objs.contour.Contours
constructor must be a dict or
an instance of :class:`plotly.graph_objs.contour.Contours`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("coloring", arg, coloring)
        self._set_property("end", arg, end)
        self._set_property("labelfont", arg, labelfont)
        self._set_property("labelformat", arg, labelformat)
        self._set_property("operation", arg, operation)
        self._set_property("showlabels", arg, showlabels)
        self._set_property("showlines", arg, showlines)
        self._set_property("size", arg, size)
        self._set_property("start", arg, start)
        self._set_property("type", arg, type)
        self._set_property("value", arg, value)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
