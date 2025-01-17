import _plotly_utils.basevalidators


class ShapesValidator(_plotly_utils.basevalidators.CompoundArrayValidator):
    def __init__(self, plotly_name="shapes", parent_name="layout", **kwargs):
        super(ShapesValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Shape"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            editable
                Determines whether the shape could be activated
                for edit or not. Has no effect when the older
                editable shapes mode is enabled via
                `config.editable` or
                `config.edits.shapePosition`.
            fillcolor
                Sets the color filling the shape's interior.
                Only applies to closed shapes.
            fillrule
                Determines which regions of complex paths
                constitute the interior. For more info please
                visit https://developer.mozilla.org/en-
                US/docs/Web/SVG/Attribute/fill-rule
            label
                :class:`plotly.graph_objects.layout.shape.Label
                ` instance or dict with compatible properties
            layer
                Specifies whether shapes are drawn below
                gridlines ("below"), between gridlines and
                traces ("between") or above traces ("above").
            legend
                Sets the reference to a legend to show this
                shape in. References to these legends are
                "legend", "legend2", "legend3", etc. Settings
                for these legends are set in the layout, under
                `layout.legend`, `layout.legend2`, etc.
            legendgroup
                Sets the legend group for this shape. Traces
                and shapes part of the same legend group
                hide/show at the same time when toggling legend
                items.
            legendgrouptitle
                :class:`plotly.graph_objects.layout.shape.Legen
                dgrouptitle` instance or dict with compatible
                properties
            legendrank
                Sets the legend rank for this shape. Items and
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
                legend for this shape.
            line
                :class:`plotly.graph_objects.layout.shape.Line`
                instance or dict with compatible properties
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
                Sets the opacity of the shape.
            path
                For `type` "path" - a valid SVG path with the
                pixel values replaced by data values in
                `xsizemode`/`ysizemode` being "scaled" and
                taken unmodified as pixels relative to
                `xanchor` and `yanchor` in case of "pixel" size
                mode. There are a few restrictions / quirks
                only absolute instructions, not relative. So
                the allowed segments are: M, L, H, V, Q, C, T,
                S, and Z arcs (A) are not allowed because
                radius rx and ry are relative. In the future we
                could consider supporting relative commands,
                but we would have to decide on how to handle
                date and log axes. Note that even as is, Q and
                C Bezier paths that are smooth on linear axes
                may not be smooth on log, and vice versa. no
                chained "polybezier" commands - specify the
                segment type for each one. On category axes,
                values are numbers scaled to the serial numbers
                of categories because using the categories
                themselves there would be no way to describe
                fractional positions On data axes: because
                space and T are both normal components of path
                strings, we can't use either to separate date
                from time parts. Therefore we'll use underscore
                for this purpose: 2015-02-21_13:45:56.789
            showlegend
                Determines whether or not this shape is shown
                in the legend.
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
            type
                Specifies the shape type to be drawn. If
                "line", a line is drawn from (`x0`,`y0`) to
                (`x1`,`y1`) with respect to the axes' sizing
                mode. If "circle", a circle is drawn from
                ((`x0`+`x1`)/2, (`y0`+`y1`)/2)) with radius
                (|(`x0`+`x1`)/2 - `x0`|, |(`y0`+`y1`)/2
                -`y0`)|) with respect to the axes' sizing mode.
                If "rect", a rectangle is drawn linking
                (`x0`,`y0`), (`x1`,`y0`), (`x1`,`y1`),
                (`x0`,`y1`), (`x0`,`y0`) with respect to the
                axes' sizing mode. If "path", draw a custom SVG
                path using `path`. with respect to the axes'
                sizing mode.
            visible
                Determines whether or not this shape is
                visible. If "legendonly", the shape is not
                drawn, but can appear as a legend item
                (provided that the legend itself is visible).
            x0
                Sets the shape's starting x position. See
                `type` and `xsizemode` for more info.
            x0shift
                Shifts `x0` away from the center of the
                category when `xref` is a "category" or
                "multicategory" axis. -0.5 corresponds to the
                start of the category and 0.5 corresponds to
                the end of the category.
            x1
                Sets the shape's end x position. See `type` and
                `xsizemode` for more info.
            x1shift
                Shifts `x1` away from the center of the
                category when `xref` is a "category" or
                "multicategory" axis. -0.5 corresponds to the
                start of the category and 0.5 corresponds to
                the end of the category.
            xanchor
                Only relevant in conjunction with `xsizemode`
                set to "pixel". Specifies the anchor point on
                the x axis to which `x0`, `x1` and x
                coordinates within `path` are relative to. E.g.
                useful to attach a pixel sized shape to a
                certain data value. No effect when `xsizemode`
                not set to "pixel".
            xref
                Sets the shape's x coordinate axis. If set to a
                x axis id (e.g. "x" or "x2"), the `x` position
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
                domain of the second x axis.
            xsizemode
                Sets the shapes's sizing mode along the x axis.
                If set to "scaled", `x0`, `x1` and x
                coordinates within `path` refer to data values
                on the x axis or a fraction of the plot area's
                width (`xref` set to "paper"). If set to
                "pixel", `xanchor` specifies the x position in
                terms of data or plot fraction but `x0`, `x1`
                and x coordinates within `path` are pixels
                relative to `xanchor`. This way, the shape can
                have a fixed width while maintaining a position
                relative to data or plot fraction.
            y0
                Sets the shape's starting y position. See
                `type` and `ysizemode` for more info.
            y0shift
                Shifts `y0` away from the center of the
                category when `yref` is a "category" or
                "multicategory" axis. -0.5 corresponds to the
                start of the category and 0.5 corresponds to
                the end of the category.
            y1
                Sets the shape's end y position. See `type` and
                `ysizemode` for more info.
            y1shift
                Shifts `y1` away from the center of the
                category when `yref` is a "category" or
                "multicategory" axis. -0.5 corresponds to the
                start of the category and 0.5 corresponds to
                the end of the category.
            yanchor
                Only relevant in conjunction with `ysizemode`
                set to "pixel". Specifies the anchor point on
                the y axis to which `y0`, `y1` and y
                coordinates within `path` are relative to. E.g.
                useful to attach a pixel sized shape to a
                certain data value. No effect when `ysizemode`
                not set to "pixel".
            yref
                Sets the shape's y coordinate axis. If set to a
                y axis id (e.g. "y" or "y2"), the `y` position
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
                domain of the second y axis.
            ysizemode
                Sets the shapes's sizing mode along the y axis.
                If set to "scaled", `y0`, `y1` and y
                coordinates within `path` refer to data values
                on the y axis or a fraction of the plot area's
                height (`yref` set to "paper"). If set to
                "pixel", `yanchor` specifies the y position in
                terms of data or plot fraction but `y0`, `y1`
                and y coordinates within `path` are pixels
                relative to `yanchor`. This way, the shape can
                have a fixed height while maintaining a
                position relative to data or plot fraction.
""",
            ),
            **kwargs,
        )
