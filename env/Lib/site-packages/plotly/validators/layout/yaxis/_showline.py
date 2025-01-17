import _plotly_utils.basevalidators


class ShowlineValidator(_plotly_utils.basevalidators.BooleanValidator):
    def __init__(self, plotly_name="showline", parent_name="layout.yaxis", **kwargs):
        super(ShowlineValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "ticks+layoutstyle"),
            **kwargs,
        )
