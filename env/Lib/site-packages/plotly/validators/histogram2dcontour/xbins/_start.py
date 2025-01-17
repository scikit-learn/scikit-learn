import _plotly_utils.basevalidators


class StartValidator(_plotly_utils.basevalidators.AnyValidator):
    def __init__(
        self, plotly_name="start", parent_name="histogram2dcontour.xbins", **kwargs
    ):
        super(StartValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
