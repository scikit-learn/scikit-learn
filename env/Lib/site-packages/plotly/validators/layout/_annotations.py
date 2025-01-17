import _plotly_utils.basevalidators


class AnnotationsValidator(_plotly_utils.basevalidators.CompoundArrayValidator):
    def __init__(self, plotly_name="annotations", parent_name="layout", **kwargs):
        super(AnnotationsValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Annotation"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            align
                Sets the horizontal alignment of the `text`
                within the box. Has an effect only if `text`
                spans two or more lines (i.e. `text` contains
                one or more <br> HTML tags) or if an explicit
                width is set to override the text width.
            arrowcolor
                Sets the color of the annotation arrow.
            arrowhead
                Sets the end annotation arrow head style.
            arrowside
                Sets the annotation arrow head position.
            arrowsize
                Sets the size of the end annotation arrow head,
                relative to `arrowwidth`. A value of 1
                (default) gives a head about 3x as wide as the
                line.
            arrowwidth
                Sets the width (in px) of annotation arrow
                line.
            ax
                Sets the x component of the arrow tail about
                the arrow head. If `axref` is `pixel`, a
                positive (negative) component corresponds to an
                arrow pointing from right to left (left to
                right). If `axref` is not `pixel` and is
                exactly the same as `xref`, this is an absolute
                value on that axis, like `x`, specified in the
                same coordinates as `xref`.
            axref
                Indicates in what coordinates the tail of the
                annotation (ax,ay) is specified. If set to a x
                axis id (e.g. "x" or "x2"), the `x` position
                refers to a x coordinate. If set to "paper",
                the `x` position refers to the distance from
                the left of the plotting area in normalized
                coordinates where 0 (1) corresponds to the left
                (right). If set to a x axis ID followed by
                "domain" (separated by a space), the position
                behaves like for "paper", but refers to the
                distance in fractions of the domain length from
                the left of the domain of that axis: e.g., *x2
                domain* refers to the domain of the second x
                axis and a x position of 0.5 refers to the
                point between the left and the right of the
                domain of the second x axis. In order for
                absolute positioning of the arrow to work,
                "axref" must be exactly the same as "xref",
                otherwise "axref" will revert to "pixel"
                (explained next). For relative positioning,
                "axref" can be set to "pixel", in which case
                the "ax" value is specified in pixels relative
                to "x". Absolute positioning is useful for
                trendline annotations which should continue to
                indicate the correct trend when zoomed.
                Relative positioning is useful for specifying
                the text offset for an annotated point.
            ay
                Sets the y component of the arrow tail about
                the arrow head. If `ayref` is `pixel`, a
                positive (negative) component corresponds to an
                arrow pointing from bottom to top (top to
                bottom). If `ayref` is not `pixel` and is
                exactly the same as `yref`, this is an absolute
                value on that axis, like `y`, specified in the
                same coordinates as `yref`.
            ayref
                Indicates in what coordinates the tail of the
                annotation (ax,ay) is specified. If set to a y
                axis id (e.g. "y" or "y2"), the `y` position
                refers to a y coordinate. If set to "paper",
                the `y` position refers to the distance from
                the bottom of the plotting area in normalized
                coordinates where 0 (1) corresponds to the
                bottom (top). If set to a y axis ID followed by
                "domain" (separated by a space), the position
                behaves like for "paper", but refers to the
                distance in fractions of the domain length from
                the bottom of the domain of that axis: e.g.,
                *y2 domain* refers to the domain of the second
                y  axis and a y position of 0.5 refers to the
                point between the bottom and the top of the
                domain of the second y axis. In order for
                absolute positioning of the arrow to work,
                "ayref" must be exactly the same as "yref",
                otherwise "ayref" will revert to "pixel"
                (explained next). For relative positioning,
                "ayref" can be set to "pixel", in which case
                the "ay" value is specified in pixels relative
                to "y". Absolute positioning is useful for
                trendline annotations which should continue to
                indicate the correct trend when zoomed.
                Relative positioning is useful for specifying
                the text offset for an annotated point.
            bgcolor
                Sets the background color of the annotation.
            bordercolor
                Sets the color of the border enclosing the
                annotation `text`.
            borderpad
                Sets the padding (in px) between the `text` and
                the enclosing border.
            borderwidth
                Sets the width (in px) of the border enclosing
                the annotation `text`.
            captureevents
                Determines whether the annotation text box
                captures mouse move and click events, or allows
                those events to pass through to data points in
                the plot that may be behind the annotation. By
                default `captureevents` is False unless
                `hovertext` is provided. If you use the event
                `plotly_clickannotation` without `hovertext`
                you must explicitly enable `captureevents`.
            clicktoshow
                Makes this annotation respond to clicks on the
                plot. If you click a data point that exactly
                matches the `x` and `y` values of this
                annotation, and it is hidden (visible: false),
                it will appear. In "onoff" mode, you must click
                the same point again to make it disappear, so
                if you click multiple points, you can show
                multiple annotations. In "onout" mode, a click
                anywhere else in the plot (on another data
                point or not) will hide this annotation. If you
                need to show/hide this annotation in response
                to different `x` or `y` values, you can set
                `xclick` and/or `yclick`. This is useful for
                example to label the side of a bar. To label
                markers though, `standoff` is preferred over
                `xclick` and `yclick`.
            font
                Sets the annotation text font.
            height
                Sets an explicit height for the text box. null
                (default) lets the text set the box height.
                Taller text will be clipped.
            hoverlabel
                :class:`plotly.graph_objects.layout.annotation.
                Hoverlabel` instance or dict with compatible
                properties
            hovertext
                Sets text to appear when hovering over this
                annotation. If omitted or blank, no hover label
                will appear.
            name
                When used in a template, named items are
                created in the output figure in addition to any
                items the figure already has in this array. You
                can modify these items in the output figure by
                making your own item with `templateitemname`
                matching this `name` alongside your
                modifications (including `visible: false` or
                `enabled: false` to hide it). Has no effect
                outside of a template.
            opacity
                Sets the opacity of the annotation (text +
                arrow).
            showarrow
                Determines whether or not the annotation is
                drawn with an arrow. If True, `text` is placed
                near the arrow's tail. If False, `text` lines
                up with the `x` and `y` provided.
            standoff
                Sets a distance, in pixels, to move the end
                arrowhead away from the position it is pointing
                at, for example to point at the edge of a
                marker independent of zoom. Note that this
                shortens the arrow from the `ax` / `ay` vector,
                in contrast to `xshift` / `yshift` which moves
                everything by this amount.
            startarrowhead
                Sets the start annotation arrow head style.
            startarrowsize
                Sets the size of the start annotation arrow
                head, relative to `arrowwidth`. A value of 1
                (default) gives a head about 3x as wide as the
                line.
            startstandoff
                Sets a distance, in pixels, to move the start
                arrowhead away from the position it is pointing
                at, for example to point at the edge of a
                marker independent of zoom. Note that this
                shortens the arrow from the `ax` / `ay` vector,
                in contrast to `xshift` / `yshift` which moves
                everything by this amount.
            templateitemname
                Used to refer to a named item in this array in
                the template. Named items from the template
                will be created even without a matching item in
                the input figure, but you can modify one by
                making an item with `templateitemname` matching
                its `name`, alongside your modifications
                (including `visible: false` or `enabled: false`
                to hide it). If there is no template or no
                matching item, this item will be hidden unless
                you explicitly show it with `visible: true`.
            text
                Sets the text associated with this annotation.
                Plotly uses a subset of HTML tags to do things
                like newline (<br>), bold (<b></b>), italics
                (<i></i>), hyperlinks (<a href='...'></a>).
                Tags <em>, <sup>, <sub>, <s>, <u> <span> are
                also supported.
            textangle
                Sets the angle at which the `text` is drawn
                with respect to the horizontal.
            valign
                Sets the vertical alignment of the `text`
                within the box. Has an effect only if an
                explicit height is set to override the text
                height.
            visible
                Determines whether or not this annotation is
                visible.
            width
                Sets an explicit width for the text box. null
                (default) lets the text set the box width.
                Wider text will be clipped. There is no
                automatic wrapping; use <br> to start a new
                line.
            x
                Sets the annotation's x position. If the axis
                `type` is "log", then you must take the log of
                your desired range. If the axis `type` is
                "date", it should be date strings, like date
                data, though Date objects and unix milliseconds
                will be accepted and converted to strings. If
                the axis `type` is "category", it should be
                numbers, using the scale where each category is
                assigned a serial number from zero in the order
                it appears.
            xanchor
                Sets the text box's horizontal position anchor
                This anchor binds the `x` position to the
                "left", "center" or "right" of the annotation.
                For example, if `x` is set to 1, `xref` to
                "paper" and `xanchor` to "right" then the
                right-most portion of the annotation lines up
                with the right-most edge of the plotting area.
                If "auto", the anchor is equivalent to "center"
                for data-referenced annotations or if there is
                an arrow, whereas for paper-referenced with no
                arrow, the anchor picked corresponds to the
                closest side.
            xclick
                Toggle this annotation when clicking a data
                point whose `x` value is `xclick` rather than
                the annotation's `x` value.
            xref
                Sets the annotation's x coordinate axis. If set
                to a x axis id (e.g. "x" or "x2"), the `x`
                position refers to a x coordinate. If set to
                "paper", the `x` position refers to the
                distance from the left of the plotting area in
                normalized coordinates where 0 (1) corresponds
                to the left (right). If set to a x axis ID
                followed by "domain" (separated by a space),
                the position behaves like for "paper", but
                refers to the distance in fractions of the
                domain length from the left of the domain of
                that axis: e.g., *x2 domain* refers to the
                domain of the second x  axis and a x position
                of 0.5 refers to the point between the left and
                the right of the domain of the second x axis.
            xshift
                Shifts the position of the whole annotation and
                arrow to the right (positive) or left
                (negative) by this many pixels.
            y
                Sets the annotation's y position. If the axis
                `type` is "log", then you must take the log of
                your desired range. If the axis `type` is
                "date", it should be date strings, like date
                data, though Date objects and unix milliseconds
                will be accepted and converted to strings. If
                the axis `type` is "category", it should be
                numbers, using the scale where each category is
                assigned a serial number from zero in the order
                it appears.
            yanchor
                Sets the text box's vertical position anchor
                This anchor binds the `y` position to the
                "top", "middle" or "bottom" of the annotation.
                For example, if `y` is set to 1, `yref` to
                "paper" and `yanchor` to "top" then the top-
                most portion of the annotation lines up with
                the top-most edge of the plotting area. If
                "auto", the anchor is equivalent to "middle"
                for data-referenced annotations or if there is
                an arrow, whereas for paper-referenced with no
                arrow, the anchor picked corresponds to the
                closest side.
            yclick
                Toggle this annotation when clicking a data
                point whose `y` value is `yclick` rather than
                the annotation's `y` value.
            yref
                Sets the annotation's y coordinate axis. If set
                to a y axis id (e.g. "y" or "y2"), the `y`
                position refers to a y coordinate. If set to
                "paper", the `y` position refers to the
                distance from the bottom of the plotting area
                in normalized coordinates where 0 (1)
                corresponds to the bottom (top). If set to a y
                axis ID followed by "domain" (separated by a
                space), the position behaves like for "paper",
                but refers to the distance in fractions of the
                domain length from the bottom of the domain of
                that axis: e.g., *y2 domain* refers to the
                domain of the second y  axis and a y position
                of 0.5 refers to the point between the bottom
                and the top of the domain of the second y axis.
            yshift
                Shifts the position of the whole annotation and
                arrow up (positive) or down (negative) by this
                many pixels.
""",
            ),
            **kwargs,
        )
