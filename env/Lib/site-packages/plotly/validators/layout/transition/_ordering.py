import _plotly_utils.basevalidators


class OrderingValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(
        self, plotly_name="ordering", parent_name="layout.transition", **kwargs
    ):
        super(OrderingValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            values=kwargs.pop("values", ["layout first", "traces first"]),
            **kwargs,
        )
