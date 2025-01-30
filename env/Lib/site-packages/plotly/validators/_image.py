import _plotly_utils.basevalidators


class ImageValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="image", parent_name="", **kwargs):
        super(ImageValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Image"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            colormodel
                Color model used to map the numerical color
                components described in `z` into colors. If
                `source` is specified, this attribute will be
                set to `rgba256` otherwise it defaults to
                `rgb`.
            customdata
                Assigns extra data each datum. This may be
                useful when listening to hover, click and
                selection events. Note that, "scatter" traces
                also appends customdata items in the markers
                DOM elements
            customdatasrc
                Sets the source reference on Chart Studio Cloud
                for `customdata`.
            dx
                Set the pixel's horizontal size.
            dy
                Set the pixel's vertical size
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
                :class:`plotly.graph_objects.image.Hoverlabel`
                instance or dict with compatible properties
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
                variables `z`, `color` and `colormodel`.
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
            legendgrouptitle
                :class:`plotly.graph_objects.image.Legendgroupt
                itle` instance or dict with compatible
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
            source
                Specifies the data URI of the image to be
                visualized. The URI consists of
                "data:image/[<media subtype>][;base64],<data>"
            stream
                :class:`plotly.graph_objects.image.Stream`
                instance or dict with compatible properties
            text
                Sets the text elements associated with each z
                value.
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
            x0
                Set the image's x position. The left edge of
                the image (or the right edge if the x axis is
                reversed or dx is negative) will be found at
                xmin=x0-dx/2
            xaxis
                Sets a reference between this trace's x
                coordinates and a 2D cartesian x axis. If "x"
                (the default value), the x coordinates refer to
                `layout.xaxis`. If "x2", the x coordinates
                refer to `layout.xaxis2`, and so on.
            y0
                Set the image's y position. The top edge of the
                image (or the bottom edge if the y axis is NOT
                reversed or if dy is negative) will be found at
                ymin=y0-dy/2. By default when an image trace is
                included, the y axis will be reversed so that
                the image is right-side-up, but you can disable
                this by setting yaxis.autorange=true or by
                providing an explicit y axis range.
            yaxis
                Sets a reference between this trace's y
                coordinates and a 2D cartesian y axis. If "y"
                (the default value), the y coordinates refer to
                `layout.yaxis`. If "y2", the y coordinates
                refer to `layout.yaxis2`, and so on.
            z
                A 2-dimensional array in which each element is
                an array of 3 or 4 numbers representing a
                color.
            zmax
                Array defining the higher bound for each color
                component. Note that the default value will
                depend on the colormodel. For the `rgb`
                colormodel, it is [255, 255, 255]. For the
                `rgba` colormodel, it is [255, 255, 255, 1].
                For the `rgba256` colormodel, it is [255, 255,
                255, 255]. For the `hsl` colormodel, it is
                [360, 100, 100]. For the `hsla` colormodel, it
                is [360, 100, 100, 1].
            zmin
                Array defining the lower bound for each color
                component. Note that the default value will
                depend on the colormodel. For the `rgb`
                colormodel, it is [0, 0, 0]. For the `rgba`
                colormodel, it is [0, 0, 0, 0]. For the
                `rgba256` colormodel, it is [0, 0, 0, 0]. For
                the `hsl` colormodel, it is [0, 0, 0]. For the
                `hsla` colormodel, it is [0, 0, 0, 0].
            zorder
                Sets the layer on which this trace is
                displayed, relative to other SVG traces on the
                same subplot. SVG traces with higher `zorder`
                appear in front of those with lower `zorder`.
            zsmooth
                Picks a smoothing algorithm used to smooth `z`
                data. This only applies for image traces that
                use the `source` attribute.
            zsrc
                Sets the source reference on Chart Studio Cloud
                for `z`.
""",
            ),
            **kwargs,
        )
