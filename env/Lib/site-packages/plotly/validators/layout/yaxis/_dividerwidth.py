import _plotly_utils.basevalidators


class DividerwidthValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(
        self, plotly_name="dividerwidth", parent_name="layout.yaxis", **kwargs
    ):
        super(DividerwidthValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "ticks"),
            **kwargs,
        )
