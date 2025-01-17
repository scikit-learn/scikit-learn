import _plotly_utils.basevalidators


class SelectionsValidator(_plotly_utils.basevalidators.CompoundArrayValidator):
    def __init__(self, plotly_name="selections", parent_name="layout", **kwargs):
        super(SelectionsValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Selection"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            line
                :class:`plotly.graph_objects.layout.selection.L
                ine` instance or dict with compatible
                properties
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
                Sets the opacity of the selection.
            path
                For `type` "path" - a valid SVG path similar to
                `shapes.path` in data coordinates. Allowed
                segments are: M, L and Z.
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
                Specifies the selection type to be drawn. If
                "rect", a rectangle is drawn linking
                (`x0`,`y0`), (`x1`,`y0`), (`x1`,`y1`) and
                (`x0`,`y1`). If "path", draw a custom SVG path
                using `path`.
            x0
                Sets the selection's starting x position.
            x1
                Sets the selection's end x position.
            xref
                Sets the selection's x coordinate axis. If set
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
            y0
                Sets the selection's starting y position.
            y1
                Sets the selection's end y position.
            yref
                Sets the selection's x coordinate axis. If set
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
""",
            ),
            **kwargs,
        )
