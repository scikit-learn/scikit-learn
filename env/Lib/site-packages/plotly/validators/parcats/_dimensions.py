import _plotly_utils.basevalidators


class DimensionsValidator(_plotly_utils.basevalidators.CompoundArrayValidator):
    def __init__(self, plotly_name="dimensions", parent_name="parcats", **kwargs):
        super(DimensionsValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Dimension"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            categoryarray
                Sets the order in which categories in this
                dimension appear. Only has an effect if
                `categoryorder` is set to "array". Used with
                `categoryorder`.
            categoryarraysrc
                Sets the source reference on Chart Studio Cloud
                for `categoryarray`.
            categoryorder
                Specifies the ordering logic for the categories
                in the dimension. By default, plotly uses
                "trace", which specifies the order that is
                present in the data supplied. Set
                `categoryorder` to *category ascending* or
                *category descending* if order should be
                determined by the alphanumerical order of the
                category names. Set `categoryorder` to "array"
                to derive the ordering from the attribute
                `categoryarray`. If a category is not found in
                the `categoryarray` array, the sorting behavior
                for that attribute will be identical to the
                "trace" mode. The unspecified categories will
                follow the categories in `categoryarray`.
            displayindex
                The display index of dimension, from left to
                right, zero indexed, defaults to dimension
                index.
            label
                The shown name of the dimension.
            ticktext
                Sets alternative tick labels for the categories
                in this dimension. Only has an effect if
                `categoryorder` is set to "array". Should be an
                array the same length as `categoryarray` Used
                with `categoryorder`.
            ticktextsrc
                Sets the source reference on Chart Studio Cloud
                for `ticktext`.
            values
                Dimension values. `values[n]` represents the
                category value of the `n`th point in the
                dataset, therefore the `values` vector for all
                dimensions must be the same (longer vectors
                will be truncated).
            valuessrc
                Sets the source reference on Chart Studio Cloud
                for `values`.
            visible
                Shows the dimension when set to `true` (the
                default). Hides the dimension for `false`.
""",
            ),
            **kwargs,
        )
