import _plotly_utils.basevalidators


class DimensionsValidator(_plotly_utils.basevalidators.CompoundArrayValidator):
    def __init__(self, plotly_name="dimensions", parent_name="splom", **kwargs):
        super(DimensionsValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Dimension"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            axis
                :class:`plotly.graph_objects.splom.dimension.Ax
                is` instance or dict with compatible properties
            label
                Sets the label corresponding to this splom
                dimension.
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
            values
                Sets the dimension values to be plotted.
            valuessrc
                Sets the source reference on Chart Studio Cloud
                for `values`.
            visible
                Determines whether or not this dimension is
                shown on the graph. Note that even visible
                false dimension contribute to the default grid
                generate by this splom trace.
""",
            ),
            **kwargs,
        )
