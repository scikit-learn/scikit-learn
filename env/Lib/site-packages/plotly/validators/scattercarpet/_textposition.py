import _plotly_utils.basevalidators


class TextpositionValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(
        self, plotly_name="textposition", parent_name="scattercarpet", **kwargs
    ):
        super(TextpositionValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            array_ok=kwargs.pop("array_ok", True),
            edit_type=kwargs.pop("edit_type", "calc"),
            values=kwargs.pop(
                "values",
                [
                    "top left",
                    "top center",
                    "top right",
                    "middle left",
                    "middle center",
                    "middle right",
                    "bottom left",
                    "bottom center",
                    "bottom right",
                ],
            ),
            **kwargs,
        )
