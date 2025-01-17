import _plotly_utils.basevalidators


class SeparatorsValidator(_plotly_utils.basevalidators.StringValidator):
    def __init__(self, plotly_name="separators", parent_name="layout", **kwargs):
        super(SeparatorsValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            **kwargs,
        )
