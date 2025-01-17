import _plotly_utils.basevalidators


class NewshapeValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="newshape", parent_name="layout", **kwargs):
        super(NewshapeValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Newshape"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            drawdirection
                When `dragmode` is set to "drawrect",
                "drawline" or "drawcircle" this limits the drag
                to be horizontal, vertical or diagonal. Using
                "diagonal" there is no limit e.g. in drawing
                lines in any direction. "ortho" limits the draw
                to be either horizontal or vertical.
                "horizontal" allows horizontal extend.
                "vertical" allows vertical extend.
            fillcolor
                Sets the color filling new shapes' interior.
                Please note that if using a fillcolor with
                alpha greater than half, drag inside the active
                shape starts moving the shape underneath,
                otherwise a new shape could be started over.
            fillrule
                Determines the path's interior. For more info
                please visit https://developer.mozilla.org/en-
                US/docs/Web/SVG/Attribute/fill-rule
            label
                :class:`plotly.graph_objects.layout.newshape.La
                bel` instance or dict with compatible
                properties
            layer
                Specifies whether new shapes are drawn below
                gridlines ("below"), between gridlines and
                traces ("between") or above traces ("above").
            legend
                Sets the reference to a legend to show new
                shape in. References to these legends are
                "legend", "legend2", "legend3", etc. Settings
                for these legends are set in the layout, under
                `layout.legend`, `layout.legend2`, etc.
            legendgroup
                Sets the legend group for new shape. Traces and
                shapes part of the same legend group hide/show
                at the same time when toggling legend items.
            legendgrouptitle
                :class:`plotly.graph_objects.layout.newshape.Le
                gendgrouptitle` instance or dict with
                compatible properties
            legendrank
                Sets the legend rank for new shape. Items and
                groups with smaller ranks are presented on
                top/left side while with "reversed"
                `legend.traceorder` they are on bottom/right
                side. The default legendrank is 1000, so that
                you can use ranks less than 1000 to place
                certain items before all unranked items, and
                ranks greater than 1000 to go after all
                unranked items.
            legendwidth
                Sets the width (in px or fraction) of the
                legend for new shape.
            line
                :class:`plotly.graph_objects.layout.newshape.Li
                ne` instance or dict with compatible properties
            name
                Sets new shape name. The name appears as the
                legend item.
            opacity
                Sets the opacity of new shapes.
            showlegend
                Determines whether or not new shape is shown in
                the legend.
            visible
                Determines whether or not new shape is visible.
                If "legendonly", the shape is not drawn, but
                can appear as a legend item (provided that the
                legend itself is visible).
""",
            ),
            **kwargs,
        )
