import _plotly_utils.basevalidators


class RollValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(
        self, plotly_name="roll", parent_name="layout.geo.projection.rotation", **kwargs
    ):
        super(RollValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            **kwargs,
        )
