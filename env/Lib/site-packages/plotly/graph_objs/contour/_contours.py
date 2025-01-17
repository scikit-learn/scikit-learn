from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Contours(_BaseTraceHierarchyType):

    # class properties
    # --------------------
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

    # coloring
    # --------
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

    # end
    # ---
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

    # labelfont
    # ---------
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
        plotly.graph_objs.contour.contours.Labelfont
        """
        return self["labelfont"]

    @labelfont.setter
    def labelfont(self, val):
        self["labelfont"] = val

    # labelformat
    # -----------
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

    # operation
    # ---------
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

    # showlabels
    # ----------
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

    # showlines
    # ---------
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

    # size
    # ----
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

    # start
    # -----
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

    # type
    # ----
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

    # value
    # -----
    @property
    def value(self):
        """
        Sets the value or values of the constraint boundary. When
        `operation` is set to one of the comparison values
        (=,<,>=,>,<=) "value" is expected to be a number. When
        `operation` is set to one of the interval values
        ([],(),[),(],][,)(,](,)[) "value" is expected to be an array of
        two numbers where the first is the lower bound and the second
        is the upper bound.

        The 'value' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["value"]

    @value.setter
    def value(self, val):
        self["value"] = val

    # Self properties description
    # ---------------------------
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
            (=,<,>=,>,<=) "value" is expected to be a number. When
            `operation` is set to one of the interval values
            ([],(),[),(],][,)(,](,)[) "value" is expected to be an
            array of two numbers where the first is the lower bound
            and the second is the upper bound.
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
            (=,<,>=,>,<=) "value" is expected to be a number. When
            `operation` is set to one of the interval values
            ([],(),[),(],][,)(,](,)[) "value" is expected to be an
            array of two numbers where the first is the lower bound
            and the second is the upper bound.

        Returns
        -------
        Contours
        """
        super(Contours, self).__init__("contours")

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
The first argument to the plotly.graph_objs.contour.Contours
constructor must be a dict or
an instance of :class:`plotly.graph_objs.contour.Contours`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("coloring", None)
        _v = coloring if coloring is not None else _v
        if _v is not None:
            self["coloring"] = _v
        _v = arg.pop("end", None)
        _v = end if end is not None else _v
        if _v is not None:
            self["end"] = _v
        _v = arg.pop("labelfont", None)
        _v = labelfont if labelfont is not None else _v
        if _v is not None:
            self["labelfont"] = _v
        _v = arg.pop("labelformat", None)
        _v = labelformat if labelformat is not None else _v
        if _v is not None:
            self["labelformat"] = _v
        _v = arg.pop("operation", None)
        _v = operation if operation is not None else _v
        if _v is not None:
            self["operation"] = _v
        _v = arg.pop("showlabels", None)
        _v = showlabels if showlabels is not None else _v
        if _v is not None:
            self["showlabels"] = _v
        _v = arg.pop("showlines", None)
        _v = showlines if showlines is not None else _v
        if _v is not None:
            self["showlines"] = _v
        _v = arg.pop("size", None)
        _v = size if size is not None else _v
        if _v is not None:
            self["size"] = _v
        _v = arg.pop("start", None)
        _v = start if start is not None else _v
        if _v is not None:
            self["start"] = _v
        _v = arg.pop("type", None)
        _v = type if type is not None else _v
        if _v is not None:
            self["type"] = _v
        _v = arg.pop("value", None)
        _v = value if value is not None else _v
        if _v is not None:
            self["value"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
