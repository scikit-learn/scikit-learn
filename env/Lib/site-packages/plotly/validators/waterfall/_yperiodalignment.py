import _plotly_utils.basevalidators


class YperiodalignmentValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(
        self, plotly_name="yperiodalignment", parent_name="waterfall", **kwargs
    ):
        super(YperiodalignmentValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            values=kwargs.pop("values", ["start", "middle", "end"]),
            **kwargs,
        )
