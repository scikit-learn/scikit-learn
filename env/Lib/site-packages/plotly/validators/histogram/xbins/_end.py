import _plotly_utils.basevalidators


class EndValidator(_plotly_utils.basevalidators.AnyValidator):
    def __init__(self, plotly_name="end", parent_name="histogram.xbins", **kwargs):
        super(EndValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
