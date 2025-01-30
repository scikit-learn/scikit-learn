import _plotly_utils.basevalidators


class MaxzoomValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(
        self, plotly_name="maxzoom", parent_name="scattermap.cluster", **kwargs
    ):
        super(MaxzoomValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            max=kwargs.pop("max", 24),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
