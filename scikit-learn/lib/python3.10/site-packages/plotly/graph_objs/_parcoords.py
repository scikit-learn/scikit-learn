#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceType as _BaseTraceType
import copy as _copy


class Parcoords(_BaseTraceType):
    _parent_path_str = ""
    _path_str = "parcoords"
    _valid_props = {
        "customdata",
        "customdatasrc",
        "dimensiondefaults",
        "dimensions",
        "domain",
        "ids",
        "idssrc",
        "labelangle",
        "labelfont",
        "labelside",
        "legend",
        "legendgrouptitle",
        "legendrank",
        "legendwidth",
        "line",
        "meta",
        "metasrc",
        "name",
        "rangefont",
        "stream",
        "tickfont",
        "type",
        "uid",
        "uirevision",
        "unselected",
        "visible",
    }

    @property
    def customdata(self):
        """
        Assigns extra data each datum. This may be useful when
        listening to hover, click and selection events. Note that,
        "scatter" traces also appends customdata items in the markers
        DOM elements

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
    def dimensions(self):
        """
        The dimensions (variables) of the parallel coordinates chart.
        2..60 dimensions are supported.

        The 'dimensions' property is a tuple of instances of
        Dimension that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.parcoords.Dimension
          - A list or tuple of dicts of string/value properties that
            will be passed to the Dimension constructor

        Returns
        -------
        tuple[plotly.graph_objs.parcoords.Dimension]
        """
        return self["dimensions"]

    @dimensions.setter
    def dimensions(self, val):
        self["dimensions"] = val

    @property
    def dimensiondefaults(self):
        """
        When used in a template (as
        layout.template.data.parcoords.dimensiondefaults), sets the
        default property values to use for elements of
        parcoords.dimensions

        The 'dimensiondefaults' property is an instance of Dimension
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.parcoords.Dimension`
          - A dict of string/value properties that will be passed
            to the Dimension constructor

        Returns
        -------
        plotly.graph_objs.parcoords.Dimension
        """
        return self["dimensiondefaults"]

    @dimensiondefaults.setter
    def dimensiondefaults(self, val):
        self["dimensiondefaults"] = val

    @property
    def domain(self):
        """
        The 'domain' property is an instance of Domain
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.parcoords.Domain`
          - A dict of string/value properties that will be passed
            to the Domain constructor

        Returns
        -------
        plotly.graph_objs.parcoords.Domain
        """
        return self["domain"]

    @domain.setter
    def domain(self, val):
        self["domain"] = val

    @property
    def ids(self):
        """
        Assigns id labels to each datum. These ids for object constancy
        of data points during animation. Should be an array of strings,
        not numbers or any other type.

        The 'ids' property is an array that may be specified as a tuple,
        list, numpy array, or pandas Series

        Returns
        -------
        numpy.ndarray
        """
        return self["ids"]

    @ids.setter
    def ids(self, val):
        self["ids"] = val

    @property
    def idssrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `ids`.

        The 'idssrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["idssrc"]

    @idssrc.setter
    def idssrc(self, val):
        self["idssrc"] = val

    @property
    def labelangle(self):
        """
        Sets the angle of the labels with respect to the horizontal.
        For example, a `tickangle` of -90 draws the labels vertically.
        Tilted labels with "labelangle" may be positioned better inside
        margins when `labelposition` is set to "bottom".

        The 'labelangle' property is a angle (in degrees) that may be
        specified as a number between -180 and 180.
        Numeric values outside this range are converted to the equivalent value
        (e.g. 270 is converted to -90).

        Returns
        -------
        int|float
        """
        return self["labelangle"]

    @labelangle.setter
    def labelangle(self, val):
        self["labelangle"] = val

    @property
    def labelfont(self):
        """
        Sets the font for the `dimension` labels.

        The 'labelfont' property is an instance of Labelfont
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.parcoords.Labelfont`
          - A dict of string/value properties that will be passed
            to the Labelfont constructor

        Returns
        -------
        plotly.graph_objs.parcoords.Labelfont
        """
        return self["labelfont"]

    @labelfont.setter
    def labelfont(self, val):
        self["labelfont"] = val

    @property
    def labelside(self):
        """
        Specifies the location of the `label`. "top" positions labels
        above, next to the title "bottom" positions labels below the
        graph Tilted labels with "labelangle" may be positioned better
        inside margins when `labelposition` is set to "bottom".

        The 'labelside' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['top', 'bottom']

        Returns
        -------
        Any
        """
        return self["labelside"]

    @labelside.setter
    def labelside(self, val):
        self["labelside"] = val

    @property
    def legend(self):
        """
        Sets the reference to a legend to show this trace in.
        References to these legends are "legend", "legend2", "legend3",
        etc. Settings for these legends are set in the layout, under
        `layout.legend`, `layout.legend2`, etc.

        The 'legend' property is an identifier of a particular
        subplot, of type 'legend', that may be specified as the string 'legend'
        optionally followed by an integer >= 1
        (e.g. 'legend', 'legend1', 'legend2', 'legend3', etc.)

        Returns
        -------
        str
        """
        return self["legend"]

    @legend.setter
    def legend(self, val):
        self["legend"] = val

    @property
    def legendgrouptitle(self):
        """
        The 'legendgrouptitle' property is an instance of Legendgrouptitle
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.parcoords.Legendgrouptitle`
          - A dict of string/value properties that will be passed
            to the Legendgrouptitle constructor

        Returns
        -------
        plotly.graph_objs.parcoords.Legendgrouptitle
        """
        return self["legendgrouptitle"]

    @legendgrouptitle.setter
    def legendgrouptitle(self, val):
        self["legendgrouptitle"] = val

    @property
    def legendrank(self):
        """
        Sets the legend rank for this trace. Items and groups with
        smaller ranks are presented on top/left side while with
        "reversed" `legend.traceorder` they are on bottom/right side.
        The default legendrank is 1000, so that you can use ranks less
        than 1000 to place certain items before all unranked items, and
        ranks greater than 1000 to go after all unranked items. When
        having unranked or equal rank items shapes would be displayed
        after traces i.e. according to their order in data and layout.

        The 'legendrank' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["legendrank"]

    @legendrank.setter
    def legendrank(self, val):
        self["legendrank"] = val

    @property
    def legendwidth(self):
        """
        Sets the width (in px or fraction) of the legend for this
        trace.

        The 'legendwidth' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["legendwidth"]

    @legendwidth.setter
    def legendwidth(self, val):
        self["legendwidth"] = val

    @property
    def line(self):
        """
        The 'line' property is an instance of Line
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.parcoords.Line`
          - A dict of string/value properties that will be passed
            to the Line constructor

        Returns
        -------
        plotly.graph_objs.parcoords.Line
        """
        return self["line"]

    @line.setter
    def line(self, val):
        self["line"] = val

    @property
    def meta(self):
        """
        Assigns extra meta information associated with this trace that
        can be used in various text attributes. Attributes such as
        trace `name`, graph, axis and colorbar `title.text`, annotation
        `text` `rangeselector`, `updatemenues` and `sliders` `label`
        text all support `meta`. To access the trace `meta` values in
        an attribute in the same trace, simply use `%{meta[i]}` where
        `i` is the index or key of the `meta` item in question. To
        access trace `meta` in layout attributes, use
        `%{data[n[.meta[i]}` where `i` is the index or key of the
        `meta` and `n` is the trace index.

        The 'meta' property accepts values of any type

        Returns
        -------
        Any|numpy.ndarray
        """
        return self["meta"]

    @meta.setter
    def meta(self, val):
        self["meta"] = val

    @property
    def metasrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `meta`.

        The 'metasrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["metasrc"]

    @metasrc.setter
    def metasrc(self, val):
        self["metasrc"] = val

    @property
    def name(self):
        """
        Sets the trace name. The trace name appears as the legend item
        and on hover.

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
    def rangefont(self):
        """
        Sets the font for the `dimension` range values.

        The 'rangefont' property is an instance of Rangefont
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.parcoords.Rangefont`
          - A dict of string/value properties that will be passed
            to the Rangefont constructor

        Returns
        -------
        plotly.graph_objs.parcoords.Rangefont
        """
        return self["rangefont"]

    @rangefont.setter
    def rangefont(self, val):
        self["rangefont"] = val

    @property
    def stream(self):
        """
        The 'stream' property is an instance of Stream
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.parcoords.Stream`
          - A dict of string/value properties that will be passed
            to the Stream constructor

        Returns
        -------
        plotly.graph_objs.parcoords.Stream
        """
        return self["stream"]

    @stream.setter
    def stream(self, val):
        self["stream"] = val

    @property
    def tickfont(self):
        """
        Sets the font for the `dimension` tick values.

        The 'tickfont' property is an instance of Tickfont
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.parcoords.Tickfont`
          - A dict of string/value properties that will be passed
            to the Tickfont constructor

        Returns
        -------
        plotly.graph_objs.parcoords.Tickfont
        """
        return self["tickfont"]

    @tickfont.setter
    def tickfont(self, val):
        self["tickfont"] = val

    @property
    def uid(self):
        """
        Assign an id to this trace, Use this to provide object
        constancy between traces during animations and transitions.

        The 'uid' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["uid"]

    @uid.setter
    def uid(self, val):
        self["uid"] = val

    @property
    def uirevision(self):
        """
        Controls persistence of some user-driven changes to the trace:
        `constraintrange` in `parcoords` traces, as well as some
        `editable: true` modifications such as `name` and
        `colorbar.title`. Defaults to `layout.uirevision`. Note that
        other user-driven trace attribute changes are controlled by
        `layout` attributes: `trace.visible` is controlled by
        `layout.legend.uirevision`, `selectedpoints` is controlled by
        `layout.selectionrevision`, and `colorbar.(x|y)` (accessible
        with `config: {editable: true}`) is controlled by
        `layout.editrevision`. Trace changes are tracked by `uid`,
        which only falls back on trace index if no `uid` is provided.
        So if your app can add/remove traces before the end of the
        `data` array, such that the same trace has a different index,
        you can still preserve user-driven changes if you give each
        trace a `uid` that stays with it as it moves.

        The 'uirevision' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["uirevision"]

    @uirevision.setter
    def uirevision(self, val):
        self["uirevision"] = val

    @property
    def unselected(self):
        """
        The 'unselected' property is an instance of Unselected
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.parcoords.Unselected`
          - A dict of string/value properties that will be passed
            to the Unselected constructor

        Returns
        -------
        plotly.graph_objs.parcoords.Unselected
        """
        return self["unselected"]

    @unselected.setter
    def unselected(self, val):
        self["unselected"] = val

    @property
    def visible(self):
        """
        Determines whether or not this trace is visible. If
        "legendonly", the trace is not drawn, but can appear as a
        legend item (provided that the legend itself is visible).

        The 'visible' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                [True, False, 'legendonly']

        Returns
        -------
        Any
        """
        return self["visible"]

    @visible.setter
    def visible(self, val):
        self["visible"] = val

    @property
    def type(self):
        return self._props["type"]

    @property
    def _prop_descriptions(self):
        return """\
        customdata
            Assigns extra data each datum. This may be useful when
            listening to hover, click and selection events. Note
            that, "scatter" traces also appends customdata items in
            the markers DOM elements
        customdatasrc
            Sets the source reference on Chart Studio Cloud for
            `customdata`.
        dimensions
            The dimensions (variables) of the parallel coordinates
            chart. 2..60 dimensions are supported.
        dimensiondefaults
            When used in a template (as
            layout.template.data.parcoords.dimensiondefaults), sets
            the default property values to use for elements of
            parcoords.dimensions
        domain
            :class:`plotly.graph_objects.parcoords.Domain` instance
            or dict with compatible properties
        ids
            Assigns id labels to each datum. These ids for object
            constancy of data points during animation. Should be an
            array of strings, not numbers or any other type.
        idssrc
            Sets the source reference on Chart Studio Cloud for
            `ids`.
        labelangle
            Sets the angle of the labels with respect to the
            horizontal. For example, a `tickangle` of -90 draws the
            labels vertically. Tilted labels with "labelangle" may
            be positioned better inside margins when
            `labelposition` is set to "bottom".
        labelfont
            Sets the font for the `dimension` labels.
        labelside
            Specifies the location of the `label`. "top" positions
            labels above, next to the title "bottom" positions
            labels below the graph Tilted labels with "labelangle"
            may be positioned better inside margins when
            `labelposition` is set to "bottom".
        legend
            Sets the reference to a legend to show this trace in.
            References to these legends are "legend", "legend2",
            "legend3", etc. Settings for these legends are set in
            the layout, under `layout.legend`, `layout.legend2`,
            etc.
        legendgrouptitle
            :class:`plotly.graph_objects.parcoords.Legendgrouptitle
            ` instance or dict with compatible properties
        legendrank
            Sets the legend rank for this trace. Items and groups
            with smaller ranks are presented on top/left side while
            with "reversed" `legend.traceorder` they are on
            bottom/right side. The default legendrank is 1000, so
            that you can use ranks less than 1000 to place certain
            items before all unranked items, and ranks greater than
            1000 to go after all unranked items. When having
            unranked or equal rank items shapes would be displayed
            after traces i.e. according to their order in data and
            layout.
        legendwidth
            Sets the width (in px or fraction) of the legend for
            this trace.
        line
            :class:`plotly.graph_objects.parcoords.Line` instance
            or dict with compatible properties
        meta
            Assigns extra meta information associated with this
            trace that can be used in various text attributes.
            Attributes such as trace `name`, graph, axis and
            colorbar `title.text`, annotation `text`
            `rangeselector`, `updatemenues` and `sliders` `label`
            text all support `meta`. To access the trace `meta`
            values in an attribute in the same trace, simply use
            `%{meta[i]}` where `i` is the index or key of the
            `meta` item in question. To access trace `meta` in
            layout attributes, use `%{data[n[.meta[i]}` where `i`
            is the index or key of the `meta` and `n` is the trace
            index.
        metasrc
            Sets the source reference on Chart Studio Cloud for
            `meta`.
        name
            Sets the trace name. The trace name appears as the
            legend item and on hover.
        rangefont
            Sets the font for the `dimension` range values.
        stream
            :class:`plotly.graph_objects.parcoords.Stream` instance
            or dict with compatible properties
        tickfont
            Sets the font for the `dimension` tick values.
        uid
            Assign an id to this trace, Use this to provide object
            constancy between traces during animations and
            transitions.
        uirevision
            Controls persistence of some user-driven changes to the
            trace: `constraintrange` in `parcoords` traces, as well
            as some `editable: true` modifications such as `name`
            and `colorbar.title`. Defaults to `layout.uirevision`.
            Note that other user-driven trace attribute changes are
            controlled by `layout` attributes: `trace.visible` is
            controlled by `layout.legend.uirevision`,
            `selectedpoints` is controlled by
            `layout.selectionrevision`, and `colorbar.(x|y)`
            (accessible with `config: {editable: true}`) is
            controlled by `layout.editrevision`. Trace changes are
            tracked by `uid`, which only falls back on trace index
            if no `uid` is provided. So if your app can add/remove
            traces before the end of the `data` array, such that
            the same trace has a different index, you can still
            preserve user-driven changes if you give each trace a
            `uid` that stays with it as it moves.
        unselected
            :class:`plotly.graph_objects.parcoords.Unselected`
            instance or dict with compatible properties
        visible
            Determines whether or not this trace is visible. If
            "legendonly", the trace is not drawn, but can appear as
            a legend item (provided that the legend itself is
            visible).
        """

    def __init__(
        self,
        arg=None,
        customdata=None,
        customdatasrc=None,
        dimensions=None,
        dimensiondefaults=None,
        domain=None,
        ids=None,
        idssrc=None,
        labelangle=None,
        labelfont=None,
        labelside=None,
        legend=None,
        legendgrouptitle=None,
        legendrank=None,
        legendwidth=None,
        line=None,
        meta=None,
        metasrc=None,
        name=None,
        rangefont=None,
        stream=None,
        tickfont=None,
        uid=None,
        uirevision=None,
        unselected=None,
        visible=None,
        **kwargs,
    ):
        """
        Construct a new Parcoords object

        Parallel coordinates for multidimensional exploratory data
        analysis. The samples are specified in `dimensions`. The colors
        are set in `line.color`.

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of :class:`plotly.graph_objs.Parcoords`
        customdata
            Assigns extra data each datum. This may be useful when
            listening to hover, click and selection events. Note
            that, "scatter" traces also appends customdata items in
            the markers DOM elements
        customdatasrc
            Sets the source reference on Chart Studio Cloud for
            `customdata`.
        dimensions
            The dimensions (variables) of the parallel coordinates
            chart. 2..60 dimensions are supported.
        dimensiondefaults
            When used in a template (as
            layout.template.data.parcoords.dimensiondefaults), sets
            the default property values to use for elements of
            parcoords.dimensions
        domain
            :class:`plotly.graph_objects.parcoords.Domain` instance
            or dict with compatible properties
        ids
            Assigns id labels to each datum. These ids for object
            constancy of data points during animation. Should be an
            array of strings, not numbers or any other type.
        idssrc
            Sets the source reference on Chart Studio Cloud for
            `ids`.
        labelangle
            Sets the angle of the labels with respect to the
            horizontal. For example, a `tickangle` of -90 draws the
            labels vertically. Tilted labels with "labelangle" may
            be positioned better inside margins when
            `labelposition` is set to "bottom".
        labelfont
            Sets the font for the `dimension` labels.
        labelside
            Specifies the location of the `label`. "top" positions
            labels above, next to the title "bottom" positions
            labels below the graph Tilted labels with "labelangle"
            may be positioned better inside margins when
            `labelposition` is set to "bottom".
        legend
            Sets the reference to a legend to show this trace in.
            References to these legends are "legend", "legend2",
            "legend3", etc. Settings for these legends are set in
            the layout, under `layout.legend`, `layout.legend2`,
            etc.
        legendgrouptitle
            :class:`plotly.graph_objects.parcoords.Legendgrouptitle
            ` instance or dict with compatible properties
        legendrank
            Sets the legend rank for this trace. Items and groups
            with smaller ranks are presented on top/left side while
            with "reversed" `legend.traceorder` they are on
            bottom/right side. The default legendrank is 1000, so
            that you can use ranks less than 1000 to place certain
            items before all unranked items, and ranks greater than
            1000 to go after all unranked items. When having
            unranked or equal rank items shapes would be displayed
            after traces i.e. according to their order in data and
            layout.
        legendwidth
            Sets the width (in px or fraction) of the legend for
            this trace.
        line
            :class:`plotly.graph_objects.parcoords.Line` instance
            or dict with compatible properties
        meta
            Assigns extra meta information associated with this
            trace that can be used in various text attributes.
            Attributes such as trace `name`, graph, axis and
            colorbar `title.text`, annotation `text`
            `rangeselector`, `updatemenues` and `sliders` `label`
            text all support `meta`. To access the trace `meta`
            values in an attribute in the same trace, simply use
            `%{meta[i]}` where `i` is the index or key of the
            `meta` item in question. To access trace `meta` in
            layout attributes, use `%{data[n[.meta[i]}` where `i`
            is the index or key of the `meta` and `n` is the trace
            index.
        metasrc
            Sets the source reference on Chart Studio Cloud for
            `meta`.
        name
            Sets the trace name. The trace name appears as the
            legend item and on hover.
        rangefont
            Sets the font for the `dimension` range values.
        stream
            :class:`plotly.graph_objects.parcoords.Stream` instance
            or dict with compatible properties
        tickfont
            Sets the font for the `dimension` tick values.
        uid
            Assign an id to this trace, Use this to provide object
            constancy between traces during animations and
            transitions.
        uirevision
            Controls persistence of some user-driven changes to the
            trace: `constraintrange` in `parcoords` traces, as well
            as some `editable: true` modifications such as `name`
            and `colorbar.title`. Defaults to `layout.uirevision`.
            Note that other user-driven trace attribute changes are
            controlled by `layout` attributes: `trace.visible` is
            controlled by `layout.legend.uirevision`,
            `selectedpoints` is controlled by
            `layout.selectionrevision`, and `colorbar.(x|y)`
            (accessible with `config: {editable: true}`) is
            controlled by `layout.editrevision`. Trace changes are
            tracked by `uid`, which only falls back on trace index
            if no `uid` is provided. So if your app can add/remove
            traces before the end of the `data` array, such that
            the same trace has a different index, you can still
            preserve user-driven changes if you give each trace a
            `uid` that stays with it as it moves.
        unselected
            :class:`plotly.graph_objects.parcoords.Unselected`
            instance or dict with compatible properties
        visible
            Determines whether or not this trace is visible. If
            "legendonly", the trace is not drawn, but can appear as
            a legend item (provided that the legend itself is
            visible).

        Returns
        -------
        Parcoords
        """
        super().__init__("parcoords")
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
The first argument to the plotly.graph_objs.Parcoords
constructor must be a dict or
an instance of :class:`plotly.graph_objs.Parcoords`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("customdata", arg, customdata)
        self._set_property("customdatasrc", arg, customdatasrc)
        self._set_property("dimensions", arg, dimensions)
        self._set_property("dimensiondefaults", arg, dimensiondefaults)
        self._set_property("domain", arg, domain)
        self._set_property("ids", arg, ids)
        self._set_property("idssrc", arg, idssrc)
        self._set_property("labelangle", arg, labelangle)
        self._set_property("labelfont", arg, labelfont)
        self._set_property("labelside", arg, labelside)
        self._set_property("legend", arg, legend)
        self._set_property("legendgrouptitle", arg, legendgrouptitle)
        self._set_property("legendrank", arg, legendrank)
        self._set_property("legendwidth", arg, legendwidth)
        self._set_property("line", arg, line)
        self._set_property("meta", arg, meta)
        self._set_property("metasrc", arg, metasrc)
        self._set_property("name", arg, name)
        self._set_property("rangefont", arg, rangefont)
        self._set_property("stream", arg, stream)
        self._set_property("tickfont", arg, tickfont)
        self._set_property("uid", arg, uid)
        self._set_property("uirevision", arg, uirevision)
        self._set_property("unselected", arg, unselected)
        self._set_property("visible", arg, visible)

        self._props["type"] = "parcoords"
        arg.pop("type", None)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
