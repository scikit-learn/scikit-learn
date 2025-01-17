import _plotly_utils.basevalidators


class AutorangeValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="autorange", parent_name="layout.yaxis", **kwargs):
        super(AutorangeValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "axrange"),
            implied_edits=kwargs.pop("implied_edits", {}),
            values=kwargs.pop(
                "values",
                [True, False, "reversed", "min reversed", "max reversed", "min", "max"],
            ),
            **kwargs,
        )
