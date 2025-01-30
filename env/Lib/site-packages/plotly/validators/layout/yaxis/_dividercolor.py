import _plotly_utils.basevalidators


class DividercolorValidator(_plotly_utils.basevalidators.ColorValidator):
    def __init__(
        self, plotly_name="dividercolor", parent_name="layout.yaxis", **kwargs
    ):
        super(DividercolorValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "ticks"),
            **kwargs,
        )
