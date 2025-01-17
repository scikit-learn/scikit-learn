import _plotly_utils.basevalidators


class ShowframeValidator(_plotly_utils.basevalidators.BooleanValidator):
    def __init__(self, plotly_name="showframe", parent_name="layout.geo", **kwargs):
        super(ShowframeValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            **kwargs,
        )
