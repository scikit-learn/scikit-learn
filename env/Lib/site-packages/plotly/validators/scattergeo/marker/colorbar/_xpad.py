import _plotly_utils.basevalidators


class XpadValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(
        self, plotly_name="xpad", parent_name="scattergeo.marker.colorbar", **kwargs
    ):
        super(XpadValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
