import _plotly_utils.basevalidators


class WeightValidator(_plotly_utils.basevalidators.IntegerValidator):
    def __init__(
        self,
        plotly_name="weight",
        parent_name="scatter3d.legendgrouptitle.font",
        **kwargs,
    ):
        super(WeightValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "style"),
            extras=kwargs.pop("extras", ["normal", "bold"]),
            max=kwargs.pop("max", 1000),
            min=kwargs.pop("min", 1),
            **kwargs,
        )
