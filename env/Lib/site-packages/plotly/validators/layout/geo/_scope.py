import _plotly_utils.basevalidators


class ScopeValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="scope", parent_name="layout.geo", **kwargs):
        super(ScopeValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            values=kwargs.pop(
                "values",
                [
                    "africa",
                    "asia",
                    "europe",
                    "north america",
                    "south america",
                    "usa",
                    "world",
                ],
            ),
            **kwargs,
        )
