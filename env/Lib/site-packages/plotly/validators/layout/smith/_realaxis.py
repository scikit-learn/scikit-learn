import _plotly_utils.basevalidators


class RealaxisValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="realaxis", parent_name="layout.smith", **kwargs):
        super(RealaxisValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Realaxis"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            color
                Sets default for all colors associated with
                this axis all at once: line, font, tick, and
                grid colors. Grid color is lightened by
                blending this with the plot background
                Individual pieces can override this.
            gridcolor
                Sets the color of the grid lines.
            griddash
                Sets the dash style of lines. Set to a dash
                type string ("solid", "dot", "dash",
                "longdash", "dashdot", or "longdashdot") or a
                dash length list in px (eg "5px,10px,2px,2px").
            gridwidth
                Sets the width (in px) of the grid lines.
            hoverformat
                Sets the hover text formatting rule using d3
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
            layer
                Sets the layer on which this axis is displayed.
                If *above traces*, this axis is displayed above
                all the subplot's traces If *below traces*,
                this axis is displayed below all the subplot's
                traces, but above the grid lines. Useful when
                used together with scatter-like traces with
                `cliponaxis` set to False to show markers
                and/or text nodes above this axis.
            linecolor
                Sets the axis line color.
            linewidth
                Sets the width (in px) of the axis line.
            showgrid
                Determines whether or not grid lines are drawn.
                If True, the grid lines are drawn at every tick
                mark.
            showline
                Determines whether or not a line bounding this
                axis is drawn.
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
            side
                Determines on which side of real axis line the
                tick and tick labels appear.
            tickangle
                Sets the angle of the tick labels with respect
                to the horizontal. For example, a `tickangle`
                of -90 draws the tick labels vertically.
            tickcolor
                Sets the tick color.
            tickfont
                Sets the tick font.
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
            ticklen
                Sets the tick length (in px).
            tickprefix
                Sets a tick label prefix.
            ticks
                Determines whether ticks are drawn or not. If
                "", this axis' ticks are not drawn. If "top"
                ("bottom"), this axis' are drawn above (below)
                the axis line.
            ticksuffix
                Sets a tick label suffix.
            tickvals
                Sets the values at which ticks on this axis
                appear.
            tickvalssrc
                Sets the source reference on Chart Studio Cloud
                for `tickvals`.
            tickwidth
                Sets the tick width (in px).
            visible
                A single toggle to hide the axis while
                preserving interaction like dragging. Default
                is true when a cheater plot is present on the
                axis, otherwise false
""",
            ),
            **kwargs,
        )
