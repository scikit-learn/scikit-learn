import _plotly_utils.basevalidators


class FillmodeValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(
        self, plotly_name="fillmode", parent_name="histogram.marker.pattern", **kwargs
    ):
        super(FillmodeValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "style"),
            values=kwargs.pop("values", ["replace", "overlay"]),
            **kwargs,
        )
