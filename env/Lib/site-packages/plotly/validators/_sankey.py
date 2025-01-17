import _plotly_utils.basevalidators


class SankeyValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="sankey", parent_name="", **kwargs):
        super(SankeyValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Sankey"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            arrangement
                If value is `snap` (the default), the node
                arrangement is assisted by automatic snapping
                of elements to preserve space between nodes
                specified via `nodepad`. If value is
                `perpendicular`, the nodes can only move along
                a line perpendicular to the flow. If value is
                `freeform`, the nodes can freely move on the
                plane. If value is `fixed`, the nodes are
                stationary.
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
                :class:`plotly.graph_objects.sankey.Domain`
                instance or dict with compatible properties
            hoverinfo
                Determines which trace information appear on
                hover. If `none` or `skip` are set, no
                information is displayed upon hovering. But, if
                `none` is set, click and hover events are still
                fired. Note that this attribute is superseded
                by `node.hoverinfo` and `node.hoverinfo` for
                nodes and links respectively.
            hoverlabel
                :class:`plotly.graph_objects.sankey.Hoverlabel`
                instance or dict with compatible properties
            ids
                Assigns id labels to each datum. These ids for
                object constancy of data points during
                animation. Should be an array of strings, not
                numbers or any other type.
            idssrc
                Sets the source reference on Chart Studio Cloud
                for `ids`.
            legend
                Sets the reference to a legend to show this
                trace in. References to these legends are
                "legend", "legend2", "legend3", etc. Settings
                for these legends are set in the layout, under
                `layout.legend`, `layout.legend2`, etc.
            legendgrouptitle
                :class:`plotly.graph_objects.sankey.Legendgroup
                title` instance or dict with compatible
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
            link
                The links of the Sankey plot.
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
            node
                The nodes of the Sankey plot.
            orientation
                Sets the orientation of the Sankey diagram.
            selectedpoints
                Array containing integer indices of selected
                points. Has an effect only for traces that
                support selections. Note that an empty array
                means an empty selection where the `unselected`
                are turned on for all points, whereas, any
                other non-array values means no selection all
                where the `selected` and `unselected` styles
                have no effect.
            stream
                :class:`plotly.graph_objects.sankey.Stream`
                instance or dict with compatible properties
            textfont
                Sets the font for node labels
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
            valueformat
                Sets the value formatting rule using d3
                formatting mini-languages which are very
                similar to those in Python. For numbers, see: h
                ttps://github.com/d3/d3-format/tree/v1.4.5#d3-
                format.
            valuesuffix
                Adds a unit to follow the value in the hover
                tooltip. Add a space if a separation is
                necessary from the value.
            visible
                Determines whether or not this trace is
                visible. If "legendonly", the trace is not
                drawn, but can appear as a legend item
                (provided that the legend itself is visible).
""",
            ),
            **kwargs,
        )
