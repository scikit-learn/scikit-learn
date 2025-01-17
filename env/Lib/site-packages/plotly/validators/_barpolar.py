import _plotly_utils.basevalidators


class BarpolarValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="barpolar", parent_name="", **kwargs):
        super(BarpolarValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Barpolar"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            base
                Sets where the bar base is drawn (in radial
                axis units). In "stack" barmode, traces that
                set "base" will be excluded and drawn in
                "overlay" mode instead.
            basesrc
                Sets the source reference on Chart Studio Cloud
                for `base`.
            customdata
                Assigns extra data each datum. This may be
                useful when listening to hover, click and
                selection events. Note that, "scatter" traces
                also appends customdata items in the markers
                DOM elements
            customdatasrc
                Sets the source reference on Chart Studio Cloud
                for `customdata`.
            dr
                Sets the r coordinate step.
            dtheta
                Sets the theta coordinate step. By default, the
                `dtheta` step equals the subplot's period
                divided by the length of the `r` coordinates.
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
                :class:`plotly.graph_objects.barpolar.Hoverlabe
                l` instance or dict with compatible properties
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
                Anything contained in tag `<extra>` is
                displayed in the secondary box, for example
                "<extra>{fullData.name}</extra>". To hide the
                secondary box completely, use an empty tag
                `<extra></extra>`.
            hovertemplatesrc
                Sets the source reference on Chart Studio Cloud
                for `hovertemplate`.
            hovertext
                Same as `text`.
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
                :class:`plotly.graph_objects.barpolar.Legendgro
                uptitle` instance or dict with compatible
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
                :class:`plotly.graph_objects.barpolar.Marker`
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
            offset
                Shifts the angular position where the bar is
                drawn (in "thetatunit" units).
            offsetsrc
                Sets the source reference on Chart Studio Cloud
                for `offset`.
            opacity
                Sets the opacity of the trace.
            r
                Sets the radial coordinates
            r0
                Alternate to `r`. Builds a linear space of r
                coordinates. Use with `dr` where `r0` is the
                starting coordinate and `dr` the step.
            rsrc
                Sets the source reference on Chart Studio Cloud
                for `r`.
            selected
                :class:`plotly.graph_objects.barpolar.Selected`
                instance or dict with compatible properties
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
                :class:`plotly.graph_objects.barpolar.Stream`
                instance or dict with compatible properties
            subplot
                Sets a reference between this trace's data
                coordinates and a polar subplot. If "polar"
                (the default value), the data refer to
                `layout.polar`. If "polar2", the data refer to
                `layout.polar2`, and so on.
            text
                Sets hover text elements associated with each
                bar. If a single string, the same string
                appears over all bars. If an array of string,
                the items are mapped in order to the this
                trace's coordinates.
            textsrc
                Sets the source reference on Chart Studio Cloud
                for `text`.
            theta
                Sets the angular coordinates
            theta0
                Alternate to `theta`. Builds a linear space of
                theta coordinates. Use with `dtheta` where
                `theta0` is the starting coordinate and
                `dtheta` the step.
            thetasrc
                Sets the source reference on Chart Studio Cloud
                for `theta`.
            thetaunit
                Sets the unit of input "theta" values. Has an
                effect only when on "linear" angular axes.
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
            unselected
                :class:`plotly.graph_objects.barpolar.Unselecte
                d` instance or dict with compatible properties
            visible
                Determines whether or not this trace is
                visible. If "legendonly", the trace is not
                drawn, but can appear as a legend item
                (provided that the legend itself is visible).
            width
                Sets the bar angular width (in "thetaunit"
                units).
            widthsrc
                Sets the source reference on Chart Studio Cloud
                for `width`.
""",
            ),
            **kwargs,
        )
