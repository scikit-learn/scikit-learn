import _plotly_utils.basevalidators


class IndicatorValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="indicator", parent_name="", **kwargs):
        super(IndicatorValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Indicator"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            align
                Sets the horizontal alignment of the `text`
                within the box. Note that this attribute has no
                effect if an angular gauge is displayed: in
                this case, it is always centered
            customdata
                Assigns extra data each datum. This may be
                useful when listening to hover, click and
                selection events. Note that, "scatter" traces
                also appends customdata items in the markers
                DOM elements
            customdatasrc
                Sets the source reference on Chart Studio Cloud
                for `customdata`.
            delta
                :class:`plotly.graph_objects.indicator.Delta`
                instance or dict with compatible properties
            domain
                :class:`plotly.graph_objects.indicator.Domain`
                instance or dict with compatible properties
            gauge
                The gauge of the Indicator plot.
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
                :class:`plotly.graph_objects.indicator.Legendgr
                ouptitle` instance or dict with compatible
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
            mode
                Determines how the value is displayed on the
                graph. `number` displays the value numerically
                in text. `delta` displays the difference to a
                reference value in text. Finally, `gauge`
                displays the value graphically on an axis.
            name
                Sets the trace name. The trace name appears as
                the legend item and on hover.
            number
                :class:`plotly.graph_objects.indicator.Number`
                instance or dict with compatible properties
            stream
                :class:`plotly.graph_objects.indicator.Stream`
                instance or dict with compatible properties
            title
                :class:`plotly.graph_objects.indicator.Title`
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
            value
                Sets the number to be displayed.
            visible
                Determines whether or not this trace is
                visible. If "legendonly", the trace is not
                drawn, but can appear as a legend item
                (provided that the legend itself is visible).
""",
            ),
            **kwargs,
        )
