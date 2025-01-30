import _plotly_utils.basevalidators


class MatchesValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="matches", parent_name="layout.yaxis", **kwargs):
        super(MatchesValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            values=kwargs.pop(
                "values",
                [
                    "/^x([2-9]|[1-9][0-9]+)?( domain)?$/",
                    "/^y([2-9]|[1-9][0-9]+)?( domain)?$/",
                ],
            ),
            **kwargs,
        )
