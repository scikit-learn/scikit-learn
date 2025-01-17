import _plotly_utils.basevalidators


class ZorderValidator(_plotly_utils.basevalidators.IntegerValidator):
    def __init__(self, plotly_name="zorder", parent_name="bar", **kwargs):
        super(ZorderValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            **kwargs,
        )
