import _plotly_utils.basevalidators


class PointcloudValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="pointcloud", parent_name="", **kwargs):
        super(PointcloudValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Pointcloud"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            customdata
                Assigns extra data each datum. This may be
                useful when listening to hover, click and
                selection events. Note that, "scatter" traces
                also appends customdata items in the markers
                DOM elements
            customdatasrc
                Sets the source reference on Chart Studio Cloud
                for `customdata`.
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
                :class:`plotly.graph_objects.pointcloud.Hoverla
                bel` instance or dict with compatible
                properties
            ids
                Assigns id labels to each datum. These ids for
                object constancy of data points during
                animation. Should be an array of strings, not
                numbers or any other type.
            idssrc
                Sets the source reference on Chart Studio Cloud
                for `ids`.
            indices
                A sequential value, 0..n, supply it to avoid
                creating this array inside plotting. If
                specified, it must be a typed `Int32Array`
                array. Its length must be equal to or greater
                than the number of points. For the best
                performance and memory use, create one large
                `indices` typed array that is guaranteed to be
                at least as long as the largest number of
                points during use, and reuse it on each
                `Plotly.restyle()` call.
            indicessrc
                Sets the source reference on Chart Studio Cloud
                for `indices`.
            legend
                Sets the reference to a legend to show this
                trace in. References to these legends are
                "legend", "legend2", "legend3", etc. Settings
                for these legends are set in the layout, under
                `layout.legend`, `layout.legend2`, etc.
            legendgroup
                Sets the legend group for this trace. Traces
                and shapes part of the same legend group
                hide/show at the same time when toggling legend
                items.
            legendgrouptitle
                :class:`plotly.graph_objects.pointcloud.Legendg
                rouptitle` instance or dict with compatible
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
            marker
                :class:`plotly.graph_objects.pointcloud.Marker`
                instance or dict with compatible properties
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
            showlegend
                Determines whether or not an item corresponding
                to this trace is shown in the legend.
            stream
                :class:`plotly.graph_objects.pointcloud.Stream`
                instance or dict with compatible properties
            text
                Sets text elements associated with each (x,y)
                pair. If a single string, the same string
                appears over all the data points. If an array
                of string, the items are mapped in order to the
                this trace's (x,y) coordinates. If trace
                `hoverinfo` contains a "text" flag and
                "hovertext" is not set, these elements will be
                seen in the hover labels.
            textsrc
                Sets the source reference on Chart Studio Cloud
                for `text`.
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
                Sets the x coordinates.
            xaxis
                Sets a reference between this trace's x
                coordinates and a 2D cartesian x axis. If "x"
                (the default value), the x coordinates refer to
                `layout.xaxis`. If "x2", the x coordinates
                refer to `layout.xaxis2`, and so on.
            xbounds
                Specify `xbounds` in the shape of `[xMin, xMax]
                to avoid looping through the `xy` typed array.
                Use it in conjunction with `xy` and `ybounds`
                for the performance benefits.
            xboundssrc
                Sets the source reference on Chart Studio Cloud
                for `xbounds`.
            xsrc
                Sets the source reference on Chart Studio Cloud
                for `x`.
            xy
                Faster alternative to specifying `x` and `y`
                separately. If supplied, it must be a typed
                `Float32Array` array that represents points
                such that `xy[i * 2] = x[i]` and `xy[i * 2 + 1]
                = y[i]`
            xysrc
                Sets the source reference on Chart Studio Cloud
                for `xy`.
            y
                Sets the y coordinates.
            yaxis
                Sets a reference between this trace's y
                coordinates and a 2D cartesian y axis. If "y"
                (the default value), the y coordinates refer to
                `layout.yaxis`. If "y2", the y coordinates
                refer to `layout.yaxis2`, and so on.
            ybounds
                Specify `ybounds` in the shape of `[yMin, yMax]
                to avoid looping through the `xy` typed array.
                Use it in conjunction with `xy` and `xbounds`
                for the performance benefits.
            yboundssrc
                Sets the source reference on Chart Studio Cloud
                for `ybounds`.
            ysrc
                Sets the source reference on Chart Studio Cloud
                for `y`.
""",
            ),
            **kwargs,
        )
