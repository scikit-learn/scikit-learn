import _plotly_utils.basevalidators


class TreemapValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="treemap", parent_name="", **kwargs):
        super(TreemapValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Treemap"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            branchvalues
                Determines how the items in `values` are
                summed. When set to "total", items in `values`
                are taken to be value of all its descendants.
                When set to "remainder", items in `values`
                corresponding to the root and the branches
                sectors are taken to be the extra part not part
                of the sum of the values at their leaves.
            count
                Determines default for `values` when it is not
                provided, by inferring a 1 for each of the
                "leaves" and/or "branches", otherwise 0.
            customdata
                Assigns extra data each datum. This may be
                useful when listening to hover, click and
                selection events. Note that, "scatter" traces
                also appends customdata items in the markers
                DOM elements
            customdatasrc
                Sets the source reference on Chart Studio Cloud
                for `customdata`.
            domain
                :class:`plotly.graph_objects.treemap.Domain`
                instance or dict with compatible properties
            hoverinfo
                Determines which trace information appear on
                hover. If `none` or `skip` are set, no
                information is displayed upon hovering. But, if
                `none` is set, click and hover events are still
                fired.
            hoverinfosrc
                Sets the source reference on Chart Studio Cloud
                for `hoverinfo`.
            hoverlabel
                :class:`plotly.graph_objects.treemap.Hoverlabel
                ` instance or dict with compatible properties
            hovertemplate
                Template string used for rendering the
                information that appear on hover box. Note that
                this will override `hoverinfo`. Variables are
                inserted using %{variable}, for example "y:
                %{y}" as well as %{xother}, {%_xother},
                {%_xother_}, {%xother_}. When showing info for
                several points, "xother" will be added to those
                with different x positions from the first
                point. An underscore before or after
                "(x|y)other" will add a space on that side,
                only when this field is shown. Numbers are
                formatted using d3-format's syntax
                %{variable:d3-format}, for example "Price:
                %{y:$.2f}". https://github.com/d3/d3-
                format/tree/v1.4.5#d3-format for details on the
                formatting syntax. Dates are formatted using
                d3-time-format's syntax %{variable|d3-time-
                format}, for example "Day: %{2019-01-01|%A}".
                https://github.com/d3/d3-time-
                format/tree/v2.2.3#locale_format for details on
                the date formatting syntax. The variables
                available in `hovertemplate` are the ones
                emitted as event data described at this link
                https://plotly.com/javascript/plotlyjs-
                events/#event-data. Additionally, every
                attributes that can be specified per-point (the
                ones that are `arrayOk: true`) are available.
                Finally, the template string has access to
                variables `currentPath`, `root`, `entry`,
                `percentRoot`, `percentEntry` and
                `percentParent`. Anything contained in tag
                `<extra>` is displayed in the secondary box,
                for example "<extra>{fullData.name}</extra>".
                To hide the secondary box completely, use an
                empty tag `<extra></extra>`.
            hovertemplatesrc
                Sets the source reference on Chart Studio Cloud
                for `hovertemplate`.
            hovertext
                Sets hover text elements associated with each
                sector. If a single string, the same string
                appears for all data points. If an array of
                string, the items are mapped in order of this
                trace's sectors. To be seen, trace `hoverinfo`
                must contain a "text" flag.
            hovertextsrc
                Sets the source reference on Chart Studio Cloud
                for `hovertext`.
            ids
                Assigns id labels to each datum. These ids for
                object constancy of data points during
                animation. Should be an array of strings, not
                numbers or any other type.
            idssrc
                Sets the source reference on Chart Studio Cloud
                for `ids`.
            insidetextfont
                Sets the font used for `textinfo` lying inside
                the sector.
            labels
                Sets the labels of each of the sectors.
            labelssrc
                Sets the source reference on Chart Studio Cloud
                for `labels`.
            legend
                Sets the reference to a legend to show this
                trace in. References to these legends are
                "legend", "legend2", "legend3", etc. Settings
                for these legends are set in the layout, under
                `layout.legend`, `layout.legend2`, etc.
            legendgrouptitle
                :class:`plotly.graph_objects.treemap.Legendgrou
                ptitle` instance or dict with compatible
                properties
            legendrank
                Sets the legend rank for this trace. Items and
                groups with smaller ranks are presented on
                top/left side while with "reversed"
                `legend.traceorder` they are on bottom/right
                side. The default legendrank is 1000, so that
                you can use ranks less than 1000 to place
                certain items before all unranked items, and
                ranks greater than 1000 to go after all
                unranked items. When having unranked or equal
                rank items shapes would be displayed after
                traces i.e. according to their order in data
                and layout.
            legendwidth
                Sets the width (in px or fraction) of the
                legend for this trace.
            level
                Sets the level from which this trace hierarchy
                is rendered. Set `level` to `''` to start from
                the root node in the hierarchy. Must be an "id"
                if `ids` is filled in, otherwise plotly
                attempts to find a matching item in `labels`.
            marker
                :class:`plotly.graph_objects.treemap.Marker`
                instance or dict with compatible properties
            maxdepth
                Sets the number of rendered sectors from any
                given `level`. Set `maxdepth` to "-1" to render
                all the levels in the hierarchy.
            meta
                Assigns extra meta information associated with
                this trace that can be used in various text
                attributes. Attributes such as trace `name`,
                graph, axis and colorbar `title.text`,
                annotation `text` `rangeselector`,
                `updatemenues` and `sliders` `label` text all
                support `meta`. To access the trace `meta`
                values in an attribute in the same trace,
                simply use `%{meta[i]}` where `i` is the index
                or key of the `meta` item in question. To
                access trace `meta` in layout attributes, use
                `%{data[n[.meta[i]}` where `i` is the index or
                key of the `meta` and `n` is the trace index.
            metasrc
                Sets the source reference on Chart Studio Cloud
                for `meta`.
            name
                Sets the trace name. The trace name appears as
                the legend item and on hover.
            opacity
                Sets the opacity of the trace.
            outsidetextfont
                Sets the font used for `textinfo` lying outside
                the sector. This option refers to the root of
                the hierarchy presented on top left corner of a
                treemap graph. Please note that if a hierarchy
                has multiple root nodes, this option won't have
                any effect and `insidetextfont` would be used.
            parents
                Sets the parent sectors for each of the
                sectors. Empty string items '' are understood
                to reference the root node in the hierarchy. If
                `ids` is filled, `parents` items are understood
                to be "ids" themselves. When `ids` is not set,
                plotly attempts to find matching items in
                `labels`, but beware they must be unique.
            parentssrc
                Sets the source reference on Chart Studio Cloud
                for `parents`.
            pathbar
                :class:`plotly.graph_objects.treemap.Pathbar`
                instance or dict with compatible properties
            root
                :class:`plotly.graph_objects.treemap.Root`
                instance or dict with compatible properties
            sort
                Determines whether or not the sectors are
                reordered from largest to smallest.
            stream
                :class:`plotly.graph_objects.treemap.Stream`
                instance or dict with compatible properties
            text
                Sets text elements associated with each sector.
                If trace `textinfo` contains a "text" flag,
                these elements will be seen on the chart. If
                trace `hoverinfo` contains a "text" flag and
                "hovertext" is not set, these elements will be
                seen in the hover labels.
            textfont
                Sets the font used for `textinfo`.
            textinfo
                Determines which trace information appear on
                the graph.
            textposition
                Sets the positions of the `text` elements.
            textsrc
                Sets the source reference on Chart Studio Cloud
                for `text`.
            texttemplate
                Template string used for rendering the
                information text that appear on points. Note
                that this will override `textinfo`. Variables
                are inserted using %{variable}, for example "y:
                %{y}". Numbers are formatted using d3-format's
                syntax %{variable:d3-format}, for example
                "Price: %{y:$.2f}". https://github.com/d3/d3-
                format/tree/v1.4.5#d3-format for details on the
                formatting syntax. Dates are formatted using
                d3-time-format's syntax %{variable|d3-time-
                format}, for example "Day: %{2019-01-01|%A}".
                https://github.com/d3/d3-time-
                format/tree/v2.2.3#locale_format for details on
                the date formatting syntax. Every attributes
                that can be specified per-point (the ones that
                are `arrayOk: true`) are available. Finally,
                the template string has access to variables
                `currentPath`, `root`, `entry`, `percentRoot`,
                `percentEntry`, `percentParent`, `label` and
                `value`.
            texttemplatesrc
                Sets the source reference on Chart Studio Cloud
                for `texttemplate`.
            tiling
                :class:`plotly.graph_objects.treemap.Tiling`
                instance or dict with compatible properties
            uid
                Assign an id to this trace, Use this to provide
                object constancy between traces during
                animations and transitions.
            uirevision
                Controls persistence of some user-driven
                changes to the trace: `constraintrange` in
                `parcoords` traces, as well as some `editable:
                true` modifications such as `name` and
                `colorbar.title`. Defaults to
                `layout.uirevision`. Note that other user-
                driven trace attribute changes are controlled
                by `layout` attributes: `trace.visible` is
                controlled by `layout.legend.uirevision`,
                `selectedpoints` is controlled by
                `layout.selectionrevision`, and
                `colorbar.(x|y)` (accessible with `config:
                {editable: true}`) is controlled by
                `layout.editrevision`. Trace changes are
                tracked by `uid`, which only falls back on
                trace index if no `uid` is provided. So if your
                app can add/remove traces before the end of the
                `data` array, such that the same trace has a
                different index, you can still preserve user-
                driven changes if you give each trace a `uid`
                that stays with it as it moves.
            values
                Sets the values associated with each of the
                sectors. Use with `branchvalues` to determine
                how the values are summed.
            valuessrc
                Sets the source reference on Chart Studio Cloud
                for `values`.
            visible
                Determines whether or not this trace is
                visible. If "legendonly", the trace is not
                drawn, but can appear as a legend item
                (provided that the legend itself is visible).
""",
            ),
            **kwargs,
        )
