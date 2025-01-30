import _plotly_utils.basevalidators


class MinorgriddashValidator(_plotly_utils.basevalidators.DashValidator):
    def __init__(
        self, plotly_name="minorgriddash", parent_name="carpet.aaxis", **kwargs
    ):
        super(MinorgriddashValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            values=kwargs.pop(
                "values", ["solid", "dot", "dash", "longdash", "dashdot", "longdashdot"]
            ),
            **kwargs,
        )
