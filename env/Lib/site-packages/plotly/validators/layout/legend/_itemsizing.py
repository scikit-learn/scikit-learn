import _plotly_utils.basevalidators


class ItemsizingValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="itemsizing", parent_name="layout.legend", **kwargs):
        super(ItemsizingValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "legend"),
            values=kwargs.pop("values", ["trace", "constant"]),
            **kwargs,
        )
