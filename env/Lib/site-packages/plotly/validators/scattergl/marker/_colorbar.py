import _plotly_utils.basevalidators


class ColorbarValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(
        self, plotly_name="colorbar", parent_name="scattergl.marker", **kwargs
    ):
        super(ColorbarValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "ColorBar"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            bgcolor
                Sets the color of padded area.
            bordercolor
                Sets the axis line color.
            borderwidth
                Sets the width (in px) or the border enclosing
                this color bar.
            dtick
                Sets the step in-between ticks on this axis.
                Use with `tick0`. Must be a positive number, or
                special strings available to "log" and "date"
                axes. If the axis `type` is "log", then ticks
                are set every 10^(n*dtick) where n is the tick
                number. For example, to set a tick mark at 1,
                10, 100, 1000, ... set dtick to 1. To set tick
                marks at 1, 100, 10000, ... set dtick to 2. To
                set tick marks at 1, 5, 25, 125, 625, 3125, ...
                set dtick to log_10(5), or 0.69897000433. "log"
                has several special values; "L<f>", where `f`
                is a positive number, gives ticks linearly
                spaced in value (but not position). For example
                `tick0` = 0.1, `dtick` = "L0.5" will put ticks
                at 0.1, 0.6, 1.1, 1.6 etc. To show powers of 10
                plus small digits between, use "D1" (all
                digits) or "D2" (only 2 and 5). `tick0` is
                ignored for "D1" and "D2". If the axis `type`
                is "date", then you must convert the time to
                milliseconds. For example, to set the interval
                between ticks to one day, set `dtick` to
                86400000.0. "date" also has special values
                "M<n>" gives ticks spaced by a number of
                months. `n` must be a positive integer. To set
                ticks on the 15th of every third month, set
                `tick0` to "2000-01-15" and `dtick` to "M3". To
                set ticks every 4 years, set `dtick` to "M48"
            exponentformat
                Determines a formatting rule for the tick
                exponents. For example, consider the number
                1,000,000,000. If "none", it appears as
                1,000,000,000. If "e", 1e+9. If "E", 1E+9. If
                "power", 1x10^9 (with 9 in a super script). If
                "SI", 1G. If "B", 1B.
            labelalias
                Replacement text for specific tick or hover
                labels. For example using {US: 'USA', CA:
                'Canada'} changes US to USA and CA to Canada.
                The labels we would have shown must match the
                keys exactly, after adding any tickprefix or
                ticksuffix. For negative numbers the minus sign
                symbol used (U+2212) is wider than the regular
                ascii dash. That means you need to use −1
                instead of -1. labelalias can be used with any
                axis type, and both keys (if needed) and values
                (if desired) can include html-like tags or
                MathJax.
            len
                Sets the length of the color bar This measure
                excludes the padding of both ends. That is, the
                color bar length is this length minus the
                padding on both ends.
            lenmode
                Determines whether this color bar's length
                (i.e. the measure in the color variation
                direction) is set in units of plot "fraction"
                or in *pixels. Use `len` to set the value.
            minexponent
                Hide SI prefix for 10^n if |n| is below this
                number. This only has an effect when
                `tickformat` is "SI" or "B".
            nticks
                Specifies the maximum number of ticks for the
                particular axis. The actual number of ticks
                will be chosen automatically to be less than or
                equal to `nticks`. Has an effect only if
                `tickmode` is set to "auto".
            orientation
                Sets the orientation of the colorbar.
            outlinecolor
                Sets the axis line color.
            outlinewidth
                Sets the width (in px) of the axis line.
            separatethousands
                If "true", even 4-digit integers are separated
            showexponent
                If "all", all exponents are shown besides their
                significands. If "first", only the exponent of
                the first tick is shown. If "last", only the
                exponent of the last tick is shown. If "none",
                no exponents appear.
            showticklabels
                Determines whether or not the tick labels are
                drawn.
            showtickprefix
                If "all", all tick labels are displayed with a
                prefix. If "first", only the first tick is
                displayed with a prefix. If "last", only the
                last tick is displayed with a suffix. If
                "none", tick prefixes are hidden.
            showticksuffix
                Same as `showtickprefix` but for tick suffixes.
            thickness
                Sets the thickness of the color bar This
                measure excludes the size of the padding, ticks
                and labels.
            thicknessmode
                Determines whether this color bar's thickness
                (i.e. the measure in the constant color
                direction) is set in units of plot "fraction"
                or in "pixels". Use `thickness` to set the
                value.
            tick0
                Sets the placement of the first tick on this
                axis. Use with `dtick`. If the axis `type` is
                "log", then you must take the log of your
                starting tick (e.g. to set the starting tick to
                100, set the `tick0` to 2) except when
                `dtick`=*L<f>* (see `dtick` for more info). If
                the axis `type` is "date", it should be a date
                string, like date data. If the axis `type` is
                "category", it should be a number, using the
                scale where each category is assigned a serial
                number from zero in the order it appears.
            tickangle
                Sets the angle of the tick labels with respect
                to the horizontal. For example, a `tickangle`
                of -90 draws the tick labels vertically.
            tickcolor
                Sets the tick color.
            tickfont
                Sets the color bar's tick label font
            tickformat
                Sets the tick label formatting rule using d3
                formatting mini-languages which are very
                similar to those in Python. For numbers, see: h
                ttps://github.com/d3/d3-format/tree/v1.4.5#d3-
                format. And for dates see:
                https://github.com/d3/d3-time-
                format/tree/v2.2.3#locale_format. We add two
                items to d3's date formatter: "%h" for half of
                the year as a decimal number as well as "%{n}f"
                for fractional seconds with n digits. For
                example, *2016-10-13 09:15:23.456* with
                tickformat "%H~%M~%S.%2f" would display
                "09~15~23.46"
            tickformatstops
                A tuple of :class:`plotly.graph_objects.scatter
                gl.marker.colorbar.Tickformatstop` instances or
                dicts with compatible properties
            tickformatstopdefaults
                When used in a template (as layout.template.dat
                a.scattergl.marker.colorbar.tickformatstopdefau
                lts), sets the default property values to use
                for elements of
                scattergl.marker.colorbar.tickformatstops
            ticklabeloverflow
                Determines how we handle tick labels that would
                overflow either the graph div or the domain of
                the axis. The default value for inside tick
                labels is *hide past domain*. In other cases
                the default is *hide past div*.
            ticklabelposition
                Determines where tick labels are drawn relative
                to the ticks. Left and right options are used
                when `orientation` is "h", top and bottom when
                `orientation` is "v".
            ticklabelstep
                Sets the spacing between tick labels as
                compared to the spacing between ticks. A value
                of 1 (default) means each tick gets a label. A
                value of 2 means shows every 2nd label. A
                larger value n means only every nth tick is
                labeled. `tick0` determines which labels are
                shown. Not implemented for axes with `type`
                "log" or "multicategory", or when `tickmode` is
                "array".
            ticklen
                Sets the tick length (in px).
            tickmode
                Sets the tick mode for this axis. If "auto",
                the number of ticks is set via `nticks`. If
                "linear", the placement of the ticks is
                determined by a starting position `tick0` and a
                tick step `dtick` ("linear" is the default
                value if `tick0` and `dtick` are provided). If
                "array", the placement of the ticks is set via
                `tickvals` and the tick text is `ticktext`.
                ("array" is the default value if `tickvals` is
                provided).
            tickprefix
                Sets a tick label prefix.
            ticks
                Determines whether ticks are drawn or not. If
                "", this axis' ticks are not drawn. If
                "outside" ("inside"), this axis' are drawn
                outside (inside) the axis lines.
            ticksuffix
                Sets a tick label suffix.
            ticktext
                Sets the text displayed at the ticks position
                via `tickvals`. Only has an effect if
                `tickmode` is set to "array". Used with
                `tickvals`.
            ticktextsrc
                Sets the source reference on Chart Studio Cloud
                for `ticktext`.
            tickvals
                Sets the values at which ticks on this axis
                appear. Only has an effect if `tickmode` is set
                to "array". Used with `ticktext`.
            tickvalssrc
                Sets the source reference on Chart Studio Cloud
                for `tickvals`.
            tickwidth
                Sets the tick width (in px).
            title
                :class:`plotly.graph_objects.scattergl.marker.c
                olorbar.Title` instance or dict with compatible
                properties
            titlefont
                Deprecated: Please use
                scattergl.marker.colorbar.title.font instead.
                Sets this color bar's title font. Note that the
                title's font used to be set by the now
                deprecated `titlefont` attribute.
            titleside
                Deprecated: Please use
                scattergl.marker.colorbar.title.side instead.
                Determines the location of color bar's title
                with respect to the color bar. Defaults to
                "top" when `orientation` if "v" and  defaults
                to "right" when `orientation` if "h". Note that
                the title's location used to be set by the now
                deprecated `titleside` attribute.
            x
                Sets the x position with respect to `xref` of
                the color bar (in plot fraction). When `xref`
                is "paper", defaults to 1.02 when `orientation`
                is "v" and 0.5 when `orientation` is "h". When
                `xref` is "container", defaults to 1 when
                `orientation` is "v" and 0.5 when `orientation`
                is "h". Must be between 0 and 1 if `xref` is
                "container" and between "-2" and 3 if `xref` is
                "paper".
            xanchor
                Sets this color bar's horizontal position
                anchor. This anchor binds the `x` position to
                the "left", "center" or "right" of the color
                bar. Defaults to "left" when `orientation` is
                "v" and "center" when `orientation` is "h".
            xpad
                Sets the amount of padding (in px) along the x
                direction.
            xref
                Sets the container `x` refers to. "container"
                spans the entire `width` of the plot. "paper"
                refers to the width of the plotting area only.
            y
                Sets the y position with respect to `yref` of
                the color bar (in plot fraction). When `yref`
                is "paper", defaults to 0.5 when `orientation`
                is "v" and 1.02 when `orientation` is "h". When
                `yref` is "container", defaults to 0.5 when
                `orientation` is "v" and 1 when `orientation`
                is "h". Must be between 0 and 1 if `yref` is
                "container" and between "-2" and 3 if `yref` is
                "paper".
            yanchor
                Sets this color bar's vertical position anchor
                This anchor binds the `y` position to the
                "top", "middle" or "bottom" of the color bar.
                Defaults to "middle" when `orientation` is "v"
                and "bottom" when `orientation` is "h".
            ypad
                Sets the amount of padding (in px) along the y
                direction.
            yref
                Sets the container `y` refers to. "container"
                spans the entire `height` of the plot. "paper"
                refers to the height of the plotting area only.
""",
            ),
            **kwargs,
        )
