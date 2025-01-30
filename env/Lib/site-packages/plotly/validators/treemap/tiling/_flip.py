import _plotly_utils.basevalidators


class FlipValidator(_plotly_utils.basevalidators.FlaglistValidator):
    def __init__(self, plotly_name="flip", parent_name="treemap.tiling", **kwargs):
        super(FlipValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            flags=kwargs.pop("flags", ["x", "y"]),
            **kwargs,
        )
