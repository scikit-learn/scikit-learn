import _plotly_utils.basevalidators


class CliponaxisValidator(_plotly_utils.basevalidators.BooleanValidator):
    def __init__(
        self, plotly_name="cliponaxis", parent_name="scatterternary", **kwargs
    ):
        super(CliponaxisValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            **kwargs,
        )
