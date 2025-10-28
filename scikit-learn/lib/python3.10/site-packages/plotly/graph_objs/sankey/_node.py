#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Node(_BaseTraceHierarchyType):
    _parent_path_str = "sankey"
    _path_str = "sankey.node"
    _valid_props = {
        "align",
        "color",
        "colorsrc",
        "customdata",
        "customdatasrc",
        "groups",
        "hoverinfo",
        "hoverlabel",
        "hovertemplate",
        "hovertemplatesrc",
        "label",
        "labelsrc",
        "line",
        "pad",
        "thickness",
        "x",
        "xsrc",
        "y",
        "ysrc",
    }

    @property
    def align(self):
        """
        Sets the alignment method used to position the nodes along the
        horizontal axis.

        The 'align' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['justify', 'left', 'right', 'center']

        Returns
        -------
        Any
        """
        return self["align"]

    @align.setter
    def align(self, val):
        self["align"] = val

    @property
    def color(self):
        """
        Sets the `node` color. It can be a single value, or an array
        for specifying color for each `node`. If `node.color` is
        omitted, then the default `Plotly` color palette will be cycled
        through to have a variety of colors. These defaults are not
        fully opaque, to allow some visibility of what is beneath the
        node.

        The 'color' property is a color and may be specified as:
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
        return self["color"]

    @color.setter
    def color(self, val):
        self["color"] = val

    @property
    def colorsrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `color`.

        The 'colorsrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["colorsrc"]

    @colorsrc.setter
    def colorsrc(self, val):
        self["colorsrc"] = val

    @property
    def customdata(self):
        """
        Assigns extra data to each node.

        The 'customdata' property is an array that may be specified as a tuple,
        list, numpy array, or pandas Series

        Returns
        -------
        numpy.ndarray
        """
        return self["customdata"]

    @customdata.setter
    def customdata(self, val):
        self["customdata"] = val

    @property
    def customdatasrc(self):
        """
        Sets the source reference on Chart Studio Cloud for
        `customdata`.

        The 'customdatasrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["customdatasrc"]

    @customdatasrc.setter
    def customdatasrc(self, val):
        self["customdatasrc"] = val

    @property
    def groups(self):
        """
        Groups of nodes. Each group is defined by an array with the
        indices of the nodes it contains. Multiple groups can be
        specified.

        The 'groups' property is an info array that may be specified as:
        * a 2D list where:
          The 'groups[i][j]' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        list
        """
        return self["groups"]

    @groups.setter
    def groups(self, val):
        self["groups"] = val

    @property
    def hoverinfo(self):
        """
        Determines which trace information appear when hovering nodes.
        If `none` or `skip` are set, no information is displayed upon
        hovering. But, if `none` is set, click and hover events are
        still fired.

        The 'hoverinfo' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['all', 'none', 'skip']

        Returns
        -------
        Any
        """
        return self["hoverinfo"]

    @hoverinfo.setter
    def hoverinfo(self, val):
        self["hoverinfo"] = val

    @property
    def hoverlabel(self):
        """
        The 'hoverlabel' property is an instance of Hoverlabel
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.sankey.node.Hoverlabel`
          - A dict of string/value properties that will be passed
            to the Hoverlabel constructor

        Returns
        -------
        plotly.graph_objs.sankey.node.Hoverlabel
        """
        return self["hoverlabel"]

    @hoverlabel.setter
    def hoverlabel(self, val):
        self["hoverlabel"] = val

    @property
    def hovertemplate(self):
        """
        Template string used for rendering the information that appear
        on hover box. Note that this will override `hoverinfo`.
        Variables are inserted using %{variable}, for example "y: %{y}"
        as well as %{xother}, {%_xother}, {%_xother_}, {%xother_}. When
        showing info for several points, "xother" will be added to
        those with different x positions from the first point. An
        underscore before or after "(x|y)other" will add a space on
        that side, only when this field is shown. Numbers are formatted
        using d3-format's syntax %{variable:d3-format}, for example
        "Price: %{y:$.2f}".
        https://github.com/d3/d3-format/tree/v1.4.5#d3-format for
        details on the formatting syntax. Dates are formatted using
        d3-time-format's syntax %{variable|d3-time-format}, for example
        "Day: %{2019-01-01|%A}". https://github.com/d3/d3-time-
        format/tree/v2.2.3#locale_format for details on the date
        formatting syntax. The variables available in `hovertemplate`
        are the ones emitted as event data described at this link
        https://plotly.com/javascript/plotlyjs-events/#event-data.
        Additionally, every attributes that can be specified per-point
        (the ones that are `arrayOk: true`) are available.  Variables
        `sourceLinks` and `targetLinks` are arrays of link
        objects.Finally, the template string has access to variables
        `value` and `label`. Anything contained in tag `<extra>` is
        displayed in the secondary box, for example
        `<extra>%{fullData.name}</extra>`. To hide the secondary box
        completely, use an empty tag `<extra></extra>`.

        The 'hovertemplate' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string
          - A tuple, list, or one-dimensional numpy array of the above

        Returns
        -------
        str|numpy.ndarray
        """
        return self["hovertemplate"]

    @hovertemplate.setter
    def hovertemplate(self, val):
        self["hovertemplate"] = val

    @property
    def hovertemplatesrc(self):
        """
        Sets the source reference on Chart Studio Cloud for
        `hovertemplate`.

        The 'hovertemplatesrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["hovertemplatesrc"]

    @hovertemplatesrc.setter
    def hovertemplatesrc(self, val):
        self["hovertemplatesrc"] = val

    @property
    def label(self):
        """
        The shown name of the node.

        The 'label' property is an array that may be specified as a tuple,
        list, numpy array, or pandas Series

        Returns
        -------
        numpy.ndarray
        """
        return self["label"]

    @label.setter
    def label(self, val):
        self["label"] = val

    @property
    def labelsrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `label`.

        The 'labelsrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["labelsrc"]

    @labelsrc.setter
    def labelsrc(self, val):
        self["labelsrc"] = val

    @property
    def line(self):
        """
        The 'line' property is an instance of Line
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.sankey.node.Line`
          - A dict of string/value properties that will be passed
            to the Line constructor

        Returns
        -------
        plotly.graph_objs.sankey.node.Line
        """
        return self["line"]

    @line.setter
    def line(self, val):
        self["line"] = val

    @property
    def pad(self):
        """
        Sets the padding (in px) between the `nodes`.

        The 'pad' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["pad"]

    @pad.setter
    def pad(self, val):
        self["pad"] = val

    @property
    def thickness(self):
        """
        Sets the thickness (in px) of the `nodes`.

        The 'thickness' property is a number and may be specified as:
          - An int or float in the interval [1, inf]

        Returns
        -------
        int|float
        """
        return self["thickness"]

    @thickness.setter
    def thickness(self, val):
        self["thickness"] = val

    @property
    def x(self):
        """
        The normalized horizontal position of the node.

        The 'x' property is an array that may be specified as a tuple,
        list, numpy array, or pandas Series

        Returns
        -------
        numpy.ndarray
        """
        return self["x"]

    @x.setter
    def x(self, val):
        self["x"] = val

    @property
    def xsrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `x`.

        The 'xsrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["xsrc"]

    @xsrc.setter
    def xsrc(self, val):
        self["xsrc"] = val

    @property
    def y(self):
        """
        The normalized vertical position of the node.

        The 'y' property is an array that may be specified as a tuple,
        list, numpy array, or pandas Series

        Returns
        -------
        numpy.ndarray
        """
        return self["y"]

    @y.setter
    def y(self, val):
        self["y"] = val

    @property
    def ysrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `y`.

        The 'ysrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["ysrc"]

    @ysrc.setter
    def ysrc(self, val):
        self["ysrc"] = val

    @property
    def _prop_descriptions(self):
        return """\
        align
            Sets the alignment method used to position the nodes
            along the horizontal axis.
        color
            Sets the `node` color. It can be a single value, or an
            array for specifying color for each `node`. If
            `node.color` is omitted, then the default `Plotly`
            color palette will be cycled through to have a variety
            of colors. These defaults are not fully opaque, to
            allow some visibility of what is beneath the node.
        colorsrc
            Sets the source reference on Chart Studio Cloud for
            `color`.
        customdata
            Assigns extra data to each node.
        customdatasrc
            Sets the source reference on Chart Studio Cloud for
            `customdata`.
        groups
            Groups of nodes. Each group is defined by an array with
            the indices of the nodes it contains. Multiple groups
            can be specified.
        hoverinfo
            Determines which trace information appear when hovering
            nodes. If `none` or `skip` are set, no information is
            displayed upon hovering. But, if `none` is set, click
            and hover events are still fired.
        hoverlabel
            :class:`plotly.graph_objects.sankey.node.Hoverlabel`
            instance or dict with compatible properties
        hovertemplate
            Template string used for rendering the information that
            appear on hover box. Note that this will override
            `hoverinfo`. Variables are inserted using %{variable},
            for example "y: %{y}" as well as %{xother}, {%_xother},
            {%_xother_}, {%xother_}. When showing info for several
            points, "xother" will be added to those with different
            x positions from the first point. An underscore before
            or after "(x|y)other" will add a space on that side,
            only when this field is shown. Numbers are formatted
            using d3-format's syntax %{variable:d3-format}, for
            example "Price: %{y:$.2f}".
            https://github.com/d3/d3-format/tree/v1.4.5#d3-format
            for details on the formatting syntax. Dates are
            formatted using d3-time-format's syntax
            %{variable|d3-time-format}, for example "Day:
            %{2019-01-01|%A}". https://github.com/d3/d3-time-
            format/tree/v2.2.3#locale_format for details on the
            date formatting syntax. The variables available in
            `hovertemplate` are the ones emitted as event data
            described at this link
            https://plotly.com/javascript/plotlyjs-events/#event-
            data. Additionally, every attributes that can be
            specified per-point (the ones that are `arrayOk: true`)
            are available.  Variables `sourceLinks` and
            `targetLinks` are arrays of link objects.Finally, the
            template string has access to variables `value` and
            `label`. Anything contained in tag `<extra>` is
            displayed in the secondary box, for example
            `<extra>%{fullData.name}</extra>`. To hide the
            secondary box completely, use an empty tag
            `<extra></extra>`.
        hovertemplatesrc
            Sets the source reference on Chart Studio Cloud for
            `hovertemplate`.
        label
            The shown name of the node.
        labelsrc
            Sets the source reference on Chart Studio Cloud for
            `label`.
        line
            :class:`plotly.graph_objects.sankey.node.Line` instance
            or dict with compatible properties
        pad
            Sets the padding (in px) between the `nodes`.
        thickness
            Sets the thickness (in px) of the `nodes`.
        x
            The normalized horizontal position of the node.
        xsrc
            Sets the source reference on Chart Studio Cloud for
            `x`.
        y
            The normalized vertical position of the node.
        ysrc
            Sets the source reference on Chart Studio Cloud for
            `y`.
        """

    def __init__(
        self,
        arg=None,
        align=None,
        color=None,
        colorsrc=None,
        customdata=None,
        customdatasrc=None,
        groups=None,
        hoverinfo=None,
        hoverlabel=None,
        hovertemplate=None,
        hovertemplatesrc=None,
        label=None,
        labelsrc=None,
        line=None,
        pad=None,
        thickness=None,
        x=None,
        xsrc=None,
        y=None,
        ysrc=None,
        **kwargs,
    ):
        """
        Construct a new Node object

        The nodes of the Sankey plot.

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of :class:`plotly.graph_objs.sankey.Node`
        align
            Sets the alignment method used to position the nodes
            along the horizontal axis.
        color
            Sets the `node` color. It can be a single value, or an
            array for specifying color for each `node`. If
            `node.color` is omitted, then the default `Plotly`
            color palette will be cycled through to have a variety
            of colors. These defaults are not fully opaque, to
            allow some visibility of what is beneath the node.
        colorsrc
            Sets the source reference on Chart Studio Cloud for
            `color`.
        customdata
            Assigns extra data to each node.
        customdatasrc
            Sets the source reference on Chart Studio Cloud for
            `customdata`.
        groups
            Groups of nodes. Each group is defined by an array with
            the indices of the nodes it contains. Multiple groups
            can be specified.
        hoverinfo
            Determines which trace information appear when hovering
            nodes. If `none` or `skip` are set, no information is
            displayed upon hovering. But, if `none` is set, click
            and hover events are still fired.
        hoverlabel
            :class:`plotly.graph_objects.sankey.node.Hoverlabel`
            instance or dict with compatible properties
        hovertemplate
            Template string used for rendering the information that
            appear on hover box. Note that this will override
            `hoverinfo`. Variables are inserted using %{variable},
            for example "y: %{y}" as well as %{xother}, {%_xother},
            {%_xother_}, {%xother_}. When showing info for several
            points, "xother" will be added to those with different
            x positions from the first point. An underscore before
            or after "(x|y)other" will add a space on that side,
            only when this field is shown. Numbers are formatted
            using d3-format's syntax %{variable:d3-format}, for
            example "Price: %{y:$.2f}".
            https://github.com/d3/d3-format/tree/v1.4.5#d3-format
            for details on the formatting syntax. Dates are
            formatted using d3-time-format's syntax
            %{variable|d3-time-format}, for example "Day:
            %{2019-01-01|%A}". https://github.com/d3/d3-time-
            format/tree/v2.2.3#locale_format for details on the
            date formatting syntax. The variables available in
            `hovertemplate` are the ones emitted as event data
            described at this link
            https://plotly.com/javascript/plotlyjs-events/#event-
            data. Additionally, every attributes that can be
            specified per-point (the ones that are `arrayOk: true`)
            are available.  Variables `sourceLinks` and
            `targetLinks` are arrays of link objects.Finally, the
            template string has access to variables `value` and
            `label`. Anything contained in tag `<extra>` is
            displayed in the secondary box, for example
            `<extra>%{fullData.name}</extra>`. To hide the
            secondary box completely, use an empty tag
            `<extra></extra>`.
        hovertemplatesrc
            Sets the source reference on Chart Studio Cloud for
            `hovertemplate`.
        label
            The shown name of the node.
        labelsrc
            Sets the source reference on Chart Studio Cloud for
            `label`.
        line
            :class:`plotly.graph_objects.sankey.node.Line` instance
            or dict with compatible properties
        pad
            Sets the padding (in px) between the `nodes`.
        thickness
            Sets the thickness (in px) of the `nodes`.
        x
            The normalized horizontal position of the node.
        xsrc
            Sets the source reference on Chart Studio Cloud for
            `x`.
        y
            The normalized vertical position of the node.
        ysrc
            Sets the source reference on Chart Studio Cloud for
            `y`.

        Returns
        -------
        Node
        """
        super().__init__("node")
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
The first argument to the plotly.graph_objs.sankey.Node
constructor must be a dict or
an instance of :class:`plotly.graph_objs.sankey.Node`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("align", arg, align)
        self._set_property("color", arg, color)
        self._set_property("colorsrc", arg, colorsrc)
        self._set_property("customdata", arg, customdata)
        self._set_property("customdatasrc", arg, customdatasrc)
        self._set_property("groups", arg, groups)
        self._set_property("hoverinfo", arg, hoverinfo)
        self._set_property("hoverlabel", arg, hoverlabel)
        self._set_property("hovertemplate", arg, hovertemplate)
        self._set_property("hovertemplatesrc", arg, hovertemplatesrc)
        self._set_property("label", arg, label)
        self._set_property("labelsrc", arg, labelsrc)
        self._set_property("line", arg, line)
        self._set_property("pad", arg, pad)
        self._set_property("thickness", arg, thickness)
        self._set_property("x", arg, x)
        self._set_property("xsrc", arg, xsrc)
        self._set_property("y", arg, y)
        self._set_property("ysrc", arg, ysrc)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
