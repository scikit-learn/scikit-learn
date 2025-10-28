#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceType as _BaseTraceType
import copy as _copy


class Parcats(_BaseTraceType):
    _parent_path_str = ""
    _path_str = "parcats"
    _valid_props = {
        "arrangement",
        "bundlecolors",
        "counts",
        "countssrc",
        "dimensiondefaults",
        "dimensions",
        "domain",
        "hoverinfo",
        "hoveron",
        "hovertemplate",
        "labelfont",
        "legendgrouptitle",
        "legendwidth",
        "line",
        "meta",
        "metasrc",
        "name",
        "sortpaths",
        "stream",
        "tickfont",
        "type",
        "uid",
        "uirevision",
        "visible",
    }

    @property
    def arrangement(self):
        """
        Sets the drag interaction mode for categories and dimensions.
        If `perpendicular`, the categories can only move along a line
        perpendicular to the paths. If `freeform`, the categories can
        freely move on the plane. If `fixed`, the categories and
        dimensions are stationary.

        The 'arrangement' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['perpendicular', 'freeform', 'fixed']

        Returns
        -------
        Any
        """
        return self["arrangement"]

    @arrangement.setter
    def arrangement(self, val):
        self["arrangement"] = val

    @property
    def bundlecolors(self):
        """
        Sort paths so that like colors are bundled together within each
        category.

        The 'bundlecolors' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["bundlecolors"]

    @bundlecolors.setter
    def bundlecolors(self, val):
        self["bundlecolors"] = val

    @property
    def counts(self):
        """
        The number of observations represented by each state. Defaults
        to 1 so that each state represents one observation

        The 'counts' property is a number and may be specified as:
          - An int or float in the interval [0, inf]
          - A tuple, list, or one-dimensional numpy array of the above

        Returns
        -------
        int|float|numpy.ndarray
        """
        return self["counts"]

    @counts.setter
    def counts(self, val):
        self["counts"] = val

    @property
    def countssrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `counts`.

        The 'countssrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["countssrc"]

    @countssrc.setter
    def countssrc(self, val):
        self["countssrc"] = val

    @property
    def dimensions(self):
        """
        The dimensions (variables) of the parallel categories diagram.

        The 'dimensions' property is a tuple of instances of
        Dimension that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.parcats.Dimension
          - A list or tuple of dicts of string/value properties that
            will be passed to the Dimension constructor

        Returns
        -------
        tuple[plotly.graph_objs.parcats.Dimension]
        """
        return self["dimensions"]

    @dimensions.setter
    def dimensions(self, val):
        self["dimensions"] = val

    @property
    def dimensiondefaults(self):
        """
        When used in a template (as
        layout.template.data.parcats.dimensiondefaults), sets the
        default property values to use for elements of
        parcats.dimensions

        The 'dimensiondefaults' property is an instance of Dimension
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.parcats.Dimension`
          - A dict of string/value properties that will be passed
            to the Dimension constructor

        Returns
        -------
        plotly.graph_objs.parcats.Dimension
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
          - An instance of :class:`plotly.graph_objs.parcats.Domain`
          - A dict of string/value properties that will be passed
            to the Domain constructor

        Returns
        -------
        plotly.graph_objs.parcats.Domain
        """
        return self["domain"]

    @domain.setter
    def domain(self, val):
        self["domain"] = val

    @property
    def hoverinfo(self):
        """
        Determines which trace information appear on hover. If `none`
        or `skip` are set, no information is displayed upon hovering.
        But, if `none` is set, click and hover events are still fired.

        The 'hoverinfo' property is a flaglist and may be specified
        as a string containing:
          - Any combination of ['count', 'probability'] joined with '+' characters
            (e.g. 'count+probability')
            OR exactly one of ['all', 'none', 'skip'] (e.g. 'skip')

        Returns
        -------
        Any
        """
        return self["hoverinfo"]

    @hoverinfo.setter
    def hoverinfo(self, val):
        self["hoverinfo"] = val

    @property
    def hoveron(self):
        """
        Sets the hover interaction mode for the parcats diagram. If
        `category`, hover interaction take place per category. If
        `color`, hover interactions take place per color per category.
        If `dimension`, hover interactions take place across all
        categories per dimension.

        The 'hoveron' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['category', 'color', 'dimension']

        Returns
        -------
        Any
        """
        return self["hoveron"]

    @hoveron.setter
    def hoveron(self, val):
        self["hoveron"] = val

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
        (the ones that are `arrayOk: true`) are available.  This value
        here applies when hovering over dimensions. Note that
        "categorycount", "colorcount" and "bandcolorcount" are only
        available when `hoveron` contains the "color" flag. Finally,
        the template string has access to variables `count`,
        `probability`, `category`, `categorycount`, `colorcount` and
        `bandcolorcount`. Anything contained in tag `<extra>` is
        displayed in the secondary box, for example
        `<extra>%{fullData.name}</extra>`. To hide the secondary box
        completely, use an empty tag `<extra></extra>`.

        The 'hovertemplate' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["hovertemplate"]

    @hovertemplate.setter
    def hovertemplate(self, val):
        self["hovertemplate"] = val

    @property
    def labelfont(self):
        """
        Sets the font for the `dimension` labels.

        The 'labelfont' property is an instance of Labelfont
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.parcats.Labelfont`
          - A dict of string/value properties that will be passed
            to the Labelfont constructor

        Returns
        -------
        plotly.graph_objs.parcats.Labelfont
        """
        return self["labelfont"]

    @labelfont.setter
    def labelfont(self, val):
        self["labelfont"] = val

    @property
    def legendgrouptitle(self):
        """
        The 'legendgrouptitle' property is an instance of Legendgrouptitle
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.parcats.Legendgrouptitle`
          - A dict of string/value properties that will be passed
            to the Legendgrouptitle constructor

        Returns
        -------
        plotly.graph_objs.parcats.Legendgrouptitle
        """
        return self["legendgrouptitle"]

    @legendgrouptitle.setter
    def legendgrouptitle(self, val):
        self["legendgrouptitle"] = val

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
          - An instance of :class:`plotly.graph_objs.parcats.Line`
          - A dict of string/value properties that will be passed
            to the Line constructor

        Returns
        -------
        plotly.graph_objs.parcats.Line
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
    def sortpaths(self):
        """
        Sets the path sorting algorithm. If `forward`, sort paths based
        on dimension categories from left to right. If `backward`, sort
        paths based on dimensions categories from right to left.

        The 'sortpaths' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['forward', 'backward']

        Returns
        -------
        Any
        """
        return self["sortpaths"]

    @sortpaths.setter
    def sortpaths(self, val):
        self["sortpaths"] = val

    @property
    def stream(self):
        """
        The 'stream' property is an instance of Stream
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.parcats.Stream`
          - A dict of string/value properties that will be passed
            to the Stream constructor

        Returns
        -------
        plotly.graph_objs.parcats.Stream
        """
        return self["stream"]

    @stream.setter
    def stream(self, val):
        self["stream"] = val

    @property
    def tickfont(self):
        """
        Sets the font for the `category` labels.

        The 'tickfont' property is an instance of Tickfont
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.parcats.Tickfont`
          - A dict of string/value properties that will be passed
            to the Tickfont constructor

        Returns
        -------
        plotly.graph_objs.parcats.Tickfont
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
        arrangement
            Sets the drag interaction mode for categories and
            dimensions. If `perpendicular`, the categories can only
            move along a line perpendicular to the paths. If
            `freeform`, the categories can freely move on the
            plane. If `fixed`, the categories and dimensions are
            stationary.
        bundlecolors
            Sort paths so that like colors are bundled together
            within each category.
        counts
            The number of observations represented by each state.
            Defaults to 1 so that each state represents one
            observation
        countssrc
            Sets the source reference on Chart Studio Cloud for
            `counts`.
        dimensions
            The dimensions (variables) of the parallel categories
            diagram.
        dimensiondefaults
            When used in a template (as
            layout.template.data.parcats.dimensiondefaults), sets
            the default property values to use for elements of
            parcats.dimensions
        domain
            :class:`plotly.graph_objects.parcats.Domain` instance
            or dict with compatible properties
        hoverinfo
            Determines which trace information appear on hover. If
            `none` or `skip` are set, no information is displayed
            upon hovering. But, if `none` is set, click and hover
            events are still fired.
        hoveron
            Sets the hover interaction mode for the parcats
            diagram. If `category`, hover interaction take place
            per category. If `color`, hover interactions take place
            per color per category. If `dimension`, hover
            interactions take place across all categories per
            dimension.
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
            are available.  This value here applies when hovering
            over dimensions. Note that "categorycount",
            "colorcount" and "bandcolorcount" are only available
            when `hoveron` contains the "color" flag. Finally, the
            template string has access to variables `count`,
            `probability`, `category`, `categorycount`,
            `colorcount` and `bandcolorcount`. Anything contained
            in tag `<extra>` is displayed in the secondary box, for
            example `<extra>%{fullData.name}</extra>`. To hide the
            secondary box completely, use an empty tag
            `<extra></extra>`.
        labelfont
            Sets the font for the `dimension` labels.
        legendgrouptitle
            :class:`plotly.graph_objects.parcats.Legendgrouptitle`
            instance or dict with compatible properties
        legendwidth
            Sets the width (in px or fraction) of the legend for
            this trace.
        line
            :class:`plotly.graph_objects.parcats.Line` instance or
            dict with compatible properties
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
        sortpaths
            Sets the path sorting algorithm. If `forward`, sort
            paths based on dimension categories from left to right.
            If `backward`, sort paths based on dimensions
            categories from right to left.
        stream
            :class:`plotly.graph_objects.parcats.Stream` instance
            or dict with compatible properties
        tickfont
            Sets the font for the `category` labels.
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
        visible
            Determines whether or not this trace is visible. If
            "legendonly", the trace is not drawn, but can appear as
            a legend item (provided that the legend itself is
            visible).
        """

    def __init__(
        self,
        arg=None,
        arrangement=None,
        bundlecolors=None,
        counts=None,
        countssrc=None,
        dimensions=None,
        dimensiondefaults=None,
        domain=None,
        hoverinfo=None,
        hoveron=None,
        hovertemplate=None,
        labelfont=None,
        legendgrouptitle=None,
        legendwidth=None,
        line=None,
        meta=None,
        metasrc=None,
        name=None,
        sortpaths=None,
        stream=None,
        tickfont=None,
        uid=None,
        uirevision=None,
        visible=None,
        **kwargs,
    ):
        """
        Construct a new Parcats object

        Parallel categories diagram for multidimensional categorical
        data.

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of :class:`plotly.graph_objs.Parcats`
        arrangement
            Sets the drag interaction mode for categories and
            dimensions. If `perpendicular`, the categories can only
            move along a line perpendicular to the paths. If
            `freeform`, the categories can freely move on the
            plane. If `fixed`, the categories and dimensions are
            stationary.
        bundlecolors
            Sort paths so that like colors are bundled together
            within each category.
        counts
            The number of observations represented by each state.
            Defaults to 1 so that each state represents one
            observation
        countssrc
            Sets the source reference on Chart Studio Cloud for
            `counts`.
        dimensions
            The dimensions (variables) of the parallel categories
            diagram.
        dimensiondefaults
            When used in a template (as
            layout.template.data.parcats.dimensiondefaults), sets
            the default property values to use for elements of
            parcats.dimensions
        domain
            :class:`plotly.graph_objects.parcats.Domain` instance
            or dict with compatible properties
        hoverinfo
            Determines which trace information appear on hover. If
            `none` or `skip` are set, no information is displayed
            upon hovering. But, if `none` is set, click and hover
            events are still fired.
        hoveron
            Sets the hover interaction mode for the parcats
            diagram. If `category`, hover interaction take place
            per category. If `color`, hover interactions take place
            per color per category. If `dimension`, hover
            interactions take place across all categories per
            dimension.
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
            are available.  This value here applies when hovering
            over dimensions. Note that "categorycount",
            "colorcount" and "bandcolorcount" are only available
            when `hoveron` contains the "color" flag. Finally, the
            template string has access to variables `count`,
            `probability`, `category`, `categorycount`,
            `colorcount` and `bandcolorcount`. Anything contained
            in tag `<extra>` is displayed in the secondary box, for
            example `<extra>%{fullData.name}</extra>`. To hide the
            secondary box completely, use an empty tag
            `<extra></extra>`.
        labelfont
            Sets the font for the `dimension` labels.
        legendgrouptitle
            :class:`plotly.graph_objects.parcats.Legendgrouptitle`
            instance or dict with compatible properties
        legendwidth
            Sets the width (in px or fraction) of the legend for
            this trace.
        line
            :class:`plotly.graph_objects.parcats.Line` instance or
            dict with compatible properties
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
        sortpaths
            Sets the path sorting algorithm. If `forward`, sort
            paths based on dimension categories from left to right.
            If `backward`, sort paths based on dimensions
            categories from right to left.
        stream
            :class:`plotly.graph_objects.parcats.Stream` instance
            or dict with compatible properties
        tickfont
            Sets the font for the `category` labels.
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
        visible
            Determines whether or not this trace is visible. If
            "legendonly", the trace is not drawn, but can appear as
            a legend item (provided that the legend itself is
            visible).

        Returns
        -------
        Parcats
        """
        super().__init__("parcats")
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
The first argument to the plotly.graph_objs.Parcats
constructor must be a dict or
an instance of :class:`plotly.graph_objs.Parcats`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("arrangement", arg, arrangement)
        self._set_property("bundlecolors", arg, bundlecolors)
        self._set_property("counts", arg, counts)
        self._set_property("countssrc", arg, countssrc)
        self._set_property("dimensions", arg, dimensions)
        self._set_property("dimensiondefaults", arg, dimensiondefaults)
        self._set_property("domain", arg, domain)
        self._set_property("hoverinfo", arg, hoverinfo)
        self._set_property("hoveron", arg, hoveron)
        self._set_property("hovertemplate", arg, hovertemplate)
        self._set_property("labelfont", arg, labelfont)
        self._set_property("legendgrouptitle", arg, legendgrouptitle)
        self._set_property("legendwidth", arg, legendwidth)
        self._set_property("line", arg, line)
        self._set_property("meta", arg, meta)
        self._set_property("metasrc", arg, metasrc)
        self._set_property("name", arg, name)
        self._set_property("sortpaths", arg, sortpaths)
        self._set_property("stream", arg, stream)
        self._set_property("tickfont", arg, tickfont)
        self._set_property("uid", arg, uid)
        self._set_property("uirevision", arg, uirevision)
        self._set_property("visible", arg, visible)

        self._props["type"] = "parcats"
        arg.pop("type", None)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
