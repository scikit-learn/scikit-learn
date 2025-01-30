import _plotly_utils.basevalidators


class CarpetValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="carpet", parent_name="", **kwargs):
        super(CarpetValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Carpet"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            a
                An array containing values of the first
                parameter value
            a0
                Alternate to `a`. Builds a linear space of a
                coordinates. Use with `da` where `a0` is the
                starting coordinate and `da` the step.
            aaxis
                :class:`plotly.graph_objects.carpet.Aaxis`
                instance or dict with compatible properties
            asrc
                Sets the source reference on Chart Studio Cloud
                for `a`.
            b
                A two dimensional array of y coordinates at
                each carpet point.
            b0
                Alternate to `b`. Builds a linear space of a
                coordinates. Use with `db` where `b0` is the
                starting coordinate and `db` the step.
            baxis
                :class:`plotly.graph_objects.carpet.Baxis`
                instance or dict with compatible properties
            bsrc
                Sets the source reference on Chart Studio Cloud
                for `b`.
            carpet
                An identifier for this carpet, so that
                `scattercarpet` and `contourcarpet` traces can
                specify a carpet plot on which they lie
            cheaterslope
                The shift applied to each successive row of
                data in creating a cheater plot. Only used if
                `x` is been omitted.
            color
                Sets default for all colors associated with
                this axis all at once: line, font, tick, and
                grid colors. Grid color is lightened by
                blending this with the plot background
                Individual pieces can override this.
            customdata
                Assigns extra data each datum. This may be
                useful when listening to hover, click and
                selection events. Note that, "scatter" traces
                also appends customdata items in the markers
                DOM elements
            customdatasrc
                Sets the source reference on Chart Studio Cloud
                for `customdata`.
            da
                Sets the a coordinate step. See `a0` for more
                info.
            db
                Sets the b coordinate step. See `b0` for more
                info.
            font
                The default font used for axis & tick labels on
                this carpet
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
                :class:`plotly.graph_objects.carpet.Legendgroup
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
            stream
                :class:`plotly.graph_objects.carpet.Stream`
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
            visible
                Determines whether or not this trace is
                visible. If "legendonly", the trace is not
                drawn, but can appear as a legend item
                (provided that the legend itself is visible).
            x
                A two dimensional array of x coordinates at
                each carpet point. If omitted, the plot is a
                cheater plot and the xaxis is hidden by
                default.
            xaxis
                Sets a reference between this trace's x
                coordinates and a 2D cartesian x axis. If "x"
                (the default value), the x coordinates refer to
                `layout.xaxis`. If "x2", the x coordinates
                refer to `layout.xaxis2`, and so on.
            xsrc
                Sets the source reference on Chart Studio Cloud
                for `x`.
            y
                A two dimensional array of y coordinates at
                each carpet point.
            yaxis
                Sets a reference between this trace's y
                coordinates and a 2D cartesian y axis. If "y"
                (the default value), the y coordinates refer to
                `layout.yaxis`. If "y2", the y coordinates
                refer to `layout.yaxis2`, and so on.
            ysrc
                Sets the source reference on Chart Studio Cloud
                for `y`.
            zorder
                Sets the layer on which this trace is
                displayed, relative to other SVG traces on the
                same subplot. SVG traces with higher `zorder`
                appear in front of those with lower `zorder`.
""",
            ),
            **kwargs,
        )
