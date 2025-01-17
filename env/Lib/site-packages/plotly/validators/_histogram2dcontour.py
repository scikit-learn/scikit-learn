import _plotly_utils.basevalidators


class Histogram2DcontourValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="histogram2dcontour", parent_name="", **kwargs):
        super(Histogram2DcontourValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Histogram2dContour"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            autobinx
                Obsolete: since v1.42 each bin attribute is
                auto-determined separately and `autobinx` is
                not needed. However, we accept `autobinx: true`
                or `false` and will update `xbins` accordingly
                before deleting `autobinx` from the trace.
            autobiny
                Obsolete: since v1.42 each bin attribute is
                auto-determined separately and `autobiny` is
                not needed. However, we accept `autobiny: true`
                or `false` and will update `ybins` accordingly
                before deleting `autobiny` from the trace.
            autocolorscale
                Determines whether the colorscale is a default
                palette (`autocolorscale: true`) or the palette
                determined by `colorscale`. In case
                `colorscale` is unspecified or `autocolorscale`
                is true, the default palette will be chosen
                according to whether numbers in the `color`
                array are all positive, all negative or mixed.
            autocontour
                Determines whether or not the contour level
                attributes are picked by an algorithm. If True,
                the number of contour levels can be set in
                `ncontours`. If False, set the contour level
                attributes in `contours`.
            bingroup
                Set the `xbingroup` and `ybingroup` default
                prefix For example, setting a `bingroup` of 1
                on two histogram2d traces will make them their
                x-bins and y-bins match separately.
            coloraxis
                Sets a reference to a shared color axis.
                References to these shared color axes are
                "coloraxis", "coloraxis2", "coloraxis3", etc.
                Settings for these shared color axes are set in
                the layout, under `layout.coloraxis`,
                `layout.coloraxis2`, etc. Note that multiple
                color scales can be linked to the same color
                axis.
            colorbar
                :class:`plotly.graph_objects.histogram2dcontour
                .ColorBar` instance or dict with compatible
                properties
            colorscale
                Sets the colorscale. The colorscale must be an
                array containing arrays mapping a normalized
                value to an rgb, rgba, hex, hsl, hsv, or named
                color string. At minimum, a mapping for the
                lowest (0) and highest (1) values are required.
                For example, `[[0, 'rgb(0,0,255)'], [1,
                'rgb(255,0,0)']]`. To control the bounds of the
                colorscale in color space, use `zmin` and
                `zmax`. Alternatively, `colorscale` may be a
                palette name string of the following list: Blac
                kbody,Bluered,Blues,Cividis,Earth,Electric,Gree
                ns,Greys,Hot,Jet,Picnic,Portland,Rainbow,RdBu,R
                eds,Viridis,YlGnBu,YlOrRd.
            contours
                :class:`plotly.graph_objects.histogram2dcontour
                .Contours` instance or dict with compatible
                properties
            customdata
                Assigns extra data each datum. This may be
                useful when listening to hover, click and
                selection events. Note that, "scatter" traces
                also appends customdata items in the markers
                DOM elements
            customdatasrc
                Sets the source reference on Chart Studio Cloud
                for `customdata`.
            histfunc
                Specifies the binning function used for this
                histogram trace. If "count", the histogram
                values are computed by counting the number of
                values lying inside each bin. If "sum", "avg",
                "min", "max", the histogram values are computed
                using the sum, the average, the minimum or the
                maximum of the values lying inside each bin
                respectively.
            histnorm
                Specifies the type of normalization used for
                this histogram trace. If "", the span of each
                bar corresponds to the number of occurrences
                (i.e. the number of data points lying inside
                the bins). If "percent" / "probability", the
                span of each bar corresponds to the percentage
                / fraction of occurrences with respect to the
                total number of sample points (here, the sum of
                all bin HEIGHTS equals 100% / 1). If "density",
                the span of each bar corresponds to the number
                of occurrences in a bin divided by the size of
                the bin interval (here, the sum of all bin
                AREAS equals the total number of sample
                points). If *probability density*, the area of
                each bar corresponds to the probability that an
                event will fall into the corresponding bin
                (here, the sum of all bin AREAS equals 1).
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
                :class:`plotly.graph_objects.histogram2dcontour
                .Hoverlabel` instance or dict with compatible
                properties
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
                variable `z` Anything contained in tag
                `<extra>` is displayed in the secondary box,
                for example "<extra>{fullData.name}</extra>".
                To hide the secondary box completely, use an
                empty tag `<extra></extra>`.
            hovertemplatesrc
                Sets the source reference on Chart Studio Cloud
                for `hovertemplate`.
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
                :class:`plotly.graph_objects.histogram2dcontour
                .Legendgrouptitle` instance or dict with
                compatible properties
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
            line
                :class:`plotly.graph_objects.histogram2dcontour
                .Line` instance or dict with compatible
                properties
            marker
                :class:`plotly.graph_objects.histogram2dcontour
                .Marker` instance or dict with compatible
                properties
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
            nbinsx
                Specifies the maximum number of desired bins.
                This value will be used in an algorithm that
                will decide the optimal bin size such that the
                histogram best visualizes the distribution of
                the data. Ignored if `xbins.size` is provided.
            nbinsy
                Specifies the maximum number of desired bins.
                This value will be used in an algorithm that
                will decide the optimal bin size such that the
                histogram best visualizes the distribution of
                the data. Ignored if `ybins.size` is provided.
            ncontours
                Sets the maximum number of contour levels. The
                actual number of contours will be chosen
                automatically to be less than or equal to the
                value of `ncontours`. Has an effect only if
                `autocontour` is True or if `contours.size` is
                missing.
            opacity
                Sets the opacity of the trace.
            reversescale
                Reverses the color mapping if true. If true,
                `zmin` will correspond to the last color in the
                array and `zmax` will correspond to the first
                color.
            showlegend
                Determines whether or not an item corresponding
                to this trace is shown in the legend.
            showscale
                Determines whether or not a colorbar is
                displayed for this trace.
            stream
                :class:`plotly.graph_objects.histogram2dcontour
                .Stream` instance or dict with compatible
                properties
            textfont
                For this trace it only has an effect if
                `coloring` is set to "heatmap". Sets the text
                font.
            texttemplate
                For this trace it only has an effect if
                `coloring` is set to "heatmap". Template string
                used for rendering the information text that
                appear on points. Note that this will override
                `textinfo`. Variables are inserted using
                %{variable}, for example "y: %{y}". Numbers are
                formatted using d3-format's syntax
                %{variable:d3-format}, for example "Price:
                %{y:$.2f}". https://github.com/d3/d3-
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
                `x`, `y`, `z` and `text`.
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
                Sets the sample data to be binned on the x
                axis.
            xaxis
                Sets a reference between this trace's x
                coordinates and a 2D cartesian x axis. If "x"
                (the default value), the x coordinates refer to
                `layout.xaxis`. If "x2", the x coordinates
                refer to `layout.xaxis2`, and so on.
            xbingroup
                Set a group of histogram traces which will have
                compatible x-bin settings. Using `xbingroup`,
                histogram2d and histogram2dcontour traces  (on
                axes of the same axis type) can have compatible
                x-bin settings. Note that the same `xbingroup`
                value can be used to set (1D) histogram
                `bingroup`
            xbins
                :class:`plotly.graph_objects.histogram2dcontour
                .XBins` instance or dict with compatible
                properties
            xcalendar
                Sets the calendar system to use with `x` date
                data.
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
            xsrc
                Sets the source reference on Chart Studio Cloud
                for `x`.
            y
                Sets the sample data to be binned on the y
                axis.
            yaxis
                Sets a reference between this trace's y
                coordinates and a 2D cartesian y axis. If "y"
                (the default value), the y coordinates refer to
                `layout.yaxis`. If "y2", the y coordinates
                refer to `layout.yaxis2`, and so on.
            ybingroup
                Set a group of histogram traces which will have
                compatible y-bin settings. Using `ybingroup`,
                histogram2d and histogram2dcontour traces  (on
                axes of the same axis type) can have compatible
                y-bin settings. Note that the same `ybingroup`
                value can be used to set (1D) histogram
                `bingroup`
            ybins
                :class:`plotly.graph_objects.histogram2dcontour
                .YBins` instance or dict with compatible
                properties
            ycalendar
                Sets the calendar system to use with `y` date
                data.
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
            ysrc
                Sets the source reference on Chart Studio Cloud
                for `y`.
            z
                Sets the aggregation data.
            zauto
                Determines whether or not the color domain is
                computed with respect to the input data (here
                in `z`) or the bounds set in `zmin` and `zmax`
                Defaults to `false` when `zmin` and `zmax` are
                set by the user.
            zhoverformat
                Sets the hover text formatting rulefor `z`
                using d3 formatting mini-languages which are
                very similar to those in Python. For numbers,
                see: https://github.com/d3/d3-
                format/tree/v1.4.5#d3-format.By default the
                values are formatted using generic number
                format.
            zmax
                Sets the upper bound of the color domain. Value
                should have the same units as in `z` and if
                set, `zmin` must be set as well.
            zmid
                Sets the mid-point of the color domain by
                scaling `zmin` and/or `zmax` to be equidistant
                to this point. Value should have the same units
                as in `z`. Has no effect when `zauto` is
                `false`.
            zmin
                Sets the lower bound of the color domain. Value
                should have the same units as in `z` and if
                set, `zmax` must be set as well.
            zsrc
                Sets the source reference on Chart Studio Cloud
                for `z`.
""",
            ),
            **kwargs,
        )
