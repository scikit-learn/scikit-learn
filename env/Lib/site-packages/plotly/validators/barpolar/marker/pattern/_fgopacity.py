import _plotly_utils.basevalidators


class FgopacityValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(
        self, plotly_name="fgopacity", parent_name="barpolar.marker.pattern", **kwargs
    ):
        super(FgopacityValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "style"),
            max=kwargs.pop("max", 1),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
