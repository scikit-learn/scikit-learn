import _plotly_utils.basevalidators


class StepValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(
        self,
        plotly_name="step",
        parent_name="layout.xaxis.rangeselector.button",
        **kwargs,
    ):
        super(StepValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            values=kwargs.pop(
                "values", ["month", "year", "day", "hour", "minute", "second", "all"]
            ),
            **kwargs,
        )
