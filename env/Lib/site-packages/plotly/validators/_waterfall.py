import _plotly_utils.basevalidators


class WaterfallValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="waterfall", parent_name="", **kwargs):
        super(WaterfallValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Waterfall"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            alignmentgroup
                Set several traces linked to the same position
                axis or matching axes to the same
                alignmentgroup. This controls whether bars
                compute their positional range dependently or
                independently.
            base
                Sets where the bar base is drawn (in position
                axis units).
            cliponaxis
                Determines whether the text nodes are clipped
                about the subplot axes. To show the text nodes
                above axis lines and tick labels, make sure to
                set `xaxis.layer` and `yaxis.layer` to *below
                traces*.
            connector
                :class:`plotly.graph_objects.waterfall.Connecto
                r` instance or dict with compatible properties
            constraintext
                Constrain the size of text inside or outside a
                bar to be no larger than the bar itself.
            customdata
                Assigns extra data each datum. This may be
                useful when listening to hover, click and
                selection events. Note that, "scatter" traces
                also appends customdata items in the markers
                DOM elements
            customdatasrc
                Sets the source reference on Chart Studio Cloud
                for `customdata`.
            decreasing
                :class:`plotly.graph_objects.waterfall.Decreasi
                ng` instance or dict with compatible properties
            dx
                Sets the x coordinate step. See `x0` for more
                info.
            dy
                Sets the y coordinate step. See `y0` for more
                info.
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
                :class:`plotly.graph_objects.waterfall.Hoverlab
                el` instance or dict with compatible properties
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
                variables `initial`, `delta` and `final`.
                Anything contained in tag `<extra>` is
                displayed in the secondary box, for example
                "<extra>{fullData.name}</extra>". To hide the
                secondary box completely, use an empty tag
                `<extra></extra>`.
            hovertemplatesrc
                Sets the source reference on Chart Studio Cloud
                for `hovertemplate`.
            hovertext
                Sets hover text elements associated with each
                (x,y) pair. If a single string, the same string
                appears over all the data points. If an array
                of string, the items are mapped in order to the
                this trace's (x,y) coordinates. To be seen,
                trace `hoverinfo` must contain a "text" flag.
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
            increasing
                :class:`plotly.graph_objects.waterfall.Increasi
                ng` instance or dict with compatible properties
            insidetextanchor
                Determines if texts are kept at center or
                start/end points in `textposition` "inside"
                mode.
            insidetextfont
                Sets the font used for `text` lying inside the
                bar.
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
                :class:`plotly.graph_objects.waterfall.Legendgr
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
            measure
                An array containing types of values. By default
                the values are considered as 'relative'.
                However; it is possible to use 'total' to
                compute the sums. Also 'absolute' could be
                applied to reset the computed total or to
                declare an initial value where needed.
            measuresrc
                Sets the source reference on Chart Studio Cloud
                for `measure`.
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
            offset
                Shifts the position where the bar is drawn (in
                position axis units). In "group" barmode,
                traces that set "offset" will be excluded and
                drawn in "overlay" mode instead.
            offsetgroup
                Set several traces linked to the same position
                axis or matching axes to the same offsetgroup
                where bars of the same position coordinate will
                line up.
            offsetsrc
                Sets the source reference on Chart Studio Cloud
                for `offset`.
            opacity
                Sets the opacity of the trace.
            orientation
                Sets the orientation of the bars. With "v"
                ("h"), the value of the each bar spans along
                the vertical (horizontal).
            outsidetextfont
                Sets the font used for `text` lying outside the
                bar.
            selectedpoints
                Array containing integer indices of selected
                points. Has an effect only for traces that
                support selections. Note that an empty array
                means an empty selection where the `unselected`
                are turned on for all points, whereas, any
                other non-array values means no selection all
                where the `selected` and `unselected` styles
                have no effect.
            showlegend
                Determines whether or not an item corresponding
                to this trace is shown in the legend.
            stream
                :class:`plotly.graph_objects.waterfall.Stream`
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
            textangle
                Sets the angle of the tick labels with respect
                to the bar. For example, a `tickangle` of -90
                draws the tick labels vertically. With "auto"
                the texts may automatically be rotated to fit
                with the maximum size in bars.
            textfont
                Sets the font used for `text`.
            textinfo
                Determines which trace information appear on
                the graph. In the case of having multiple
                waterfalls, totals are computed separately (per
                trace).
            textposition
                Specifies the location of the `text`. "inside"
                positions `text` inside, next to the bar end
                (rotated and scaled if needed). "outside"
                positions `text` outside, next to the bar end
                (scaled if needed), unless there is another bar
                stacked on this one, then the text gets pushed
                inside. "auto" tries to position `text` inside
                the bar, but if the bar is too small and no bar
                is stacked on this one the text is moved
                outside. If "none", no text appears.
            textpositionsrc
                Sets the source reference on Chart Studio Cloud
                for `textposition`.
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
                `initial`, `delta`, `final` and `label`.
            texttemplatesrc
                Sets the source reference on Chart Studio Cloud
                for `texttemplate`.
            totals
                :class:`plotly.graph_objects.waterfall.Totals`
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
            width
                Sets the bar width (in position axis units).
            widthsrc
                Sets the source reference on Chart Studio Cloud
                for `width`.
            x
                Sets the x coordinates.
            x0
                Alternate to `x`. Builds a linear space of x
                coordinates. Use with `dx` where `x0` is the
                starting coordinate and `dx` the step.
            xaxis
                Sets a reference between this trace's x
                coordinates and a 2D cartesian x axis. If "x"
                (the default value), the x coordinates refer to
                `layout.xaxis`. If "x2", the x coordinates
                refer to `layout.xaxis2`, and so on.
            xhoverformat
                Sets the hover text formatting rulefor `x`
                using d3 formatting mini-languages which are
                very similar to those in Python. For numbers,
                see: https://github.com/d3/d3-
                format/tree/v1.4.5#d3-format. And for dates
                see: https://github.com/d3/d3-time-
                format/tree/v2.2.3#locale_format. We add two
                items to d3's date formatter: "%h" for half of
                the year as a decimal number as well as "%{n}f"
                for fractional seconds with n digits. For
                example, *2016-10-13 09:15:23.456* with
                tickformat "%H~%M~%S.%2f" would display
                *09~15~23.46*By default the values are
                formatted using `xaxis.hoverformat`.
            xperiod
                Only relevant when the axis `type` is "date".
                Sets the period positioning in milliseconds or
                "M<n>" on the x axis. Special values in the
                form of "M<n>" could be used to declare the
                number of months. In this case `n` must be a
                positive integer.
            xperiod0
                Only relevant when the axis `type` is "date".
                Sets the base for period positioning in
                milliseconds or date string on the x0 axis.
                When `x0period` is round number of weeks, the
                `x0period0` by default would be on a Sunday
                i.e. 2000-01-02, otherwise it would be at
                2000-01-01.
            xperiodalignment
                Only relevant when the axis `type` is "date".
                Sets the alignment of data points on the x
                axis.
            xsrc
                Sets the source reference on Chart Studio Cloud
                for `x`.
            y
                Sets the y coordinates.
            y0
                Alternate to `y`. Builds a linear space of y
                coordinates. Use with `dy` where `y0` is the
                starting coordinate and `dy` the step.
            yaxis
                Sets a reference between this trace's y
                coordinates and a 2D cartesian y axis. If "y"
                (the default value), the y coordinates refer to
                `layout.yaxis`. If "y2", the y coordinates
                refer to `layout.yaxis2`, and so on.
            yhoverformat
                Sets the hover text formatting rulefor `y`
                using d3 formatting mini-languages which are
                very similar to those in Python. For numbers,
                see: https://github.com/d3/d3-
                format/tree/v1.4.5#d3-format. And for dates
                see: https://github.com/d3/d3-time-
                format/tree/v2.2.3#locale_format. We add two
                items to d3's date formatter: "%h" for half of
                the year as a decimal number as well as "%{n}f"
                for fractional seconds with n digits. For
                example, *2016-10-13 09:15:23.456* with
                tickformat "%H~%M~%S.%2f" would display
                *09~15~23.46*By default the values are
                formatted using `yaxis.hoverformat`.
            yperiod
                Only relevant when the axis `type` is "date".
                Sets the period positioning in milliseconds or
                "M<n>" on the y axis. Special values in the
                form of "M<n>" could be used to declare the
                number of months. In this case `n` must be a
                positive integer.
            yperiod0
                Only relevant when the axis `type` is "date".
                Sets the base for period positioning in
                milliseconds or date string on the y0 axis.
                When `y0period` is round number of weeks, the
                `y0period0` by default would be on a Sunday
                i.e. 2000-01-02, otherwise it would be at
                2000-01-01.
            yperiodalignment
                Only relevant when the axis `type` is "date".
                Sets the alignment of data points on the y
                axis.
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
