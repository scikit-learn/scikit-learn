import _plotly_utils.basevalidators


class SolidityValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(
        self, plotly_name="solidity", parent_name="scatter.fillpattern", **kwargs
    ):
        super(SolidityValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            array_ok=kwargs.pop("array_ok", True),
            edit_type=kwargs.pop("edit_type", "style"),
            max=kwargs.pop("max", 1),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
