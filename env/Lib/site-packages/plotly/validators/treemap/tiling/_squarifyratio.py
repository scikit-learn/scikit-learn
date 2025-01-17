import _plotly_utils.basevalidators


class SquarifyratioValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(
        self, plotly_name="squarifyratio", parent_name="treemap.tiling", **kwargs
    ):
        super(SquarifyratioValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            min=kwargs.pop("min", 1),
            **kwargs,
        )
