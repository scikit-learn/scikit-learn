import _plotly_utils.basevalidators


class TicklabelstepValidator(_plotly_utils.basevalidators.IntegerValidator):
    def __init__(
        self,
        plotly_name="ticklabelstep",
        parent_name="treemap.marker.colorbar",
        **kwargs,
    ):
        super(TicklabelstepValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "colorbars"),
            min=kwargs.pop("min", 1),
            **kwargs,
        )
