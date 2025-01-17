import _plotly_utils.basevalidators


class PatternValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(
        self, plotly_name="pattern", parent_name="layout.xaxis.rangebreak", **kwargs
    ):
        super(PatternValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            values=kwargs.pop("values", ["day of week", "hour", ""]),
            **kwargs,
        )
