import _plotly_utils.basevalidators


class PositionValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="position", parent_name="pie.title", **kwargs):
        super(PositionValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            values=kwargs.pop(
                "values",
                [
                    "top left",
                    "top center",
                    "top right",
                    "middle center",
                    "bottom left",
                    "bottom center",
                    "bottom right",
                ],
            ),
            **kwargs,
        )
