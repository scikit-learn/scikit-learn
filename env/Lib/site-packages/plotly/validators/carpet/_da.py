import _plotly_utils.basevalidators


class DaValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="da", parent_name="carpet", **kwargs):
        super(DaValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
