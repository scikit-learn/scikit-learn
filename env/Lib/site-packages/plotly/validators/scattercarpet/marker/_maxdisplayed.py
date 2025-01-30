import _plotly_utils.basevalidators


class MaxdisplayedValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(
        self, plotly_name="maxdisplayed", parent_name="scattercarpet.marker", **kwargs
    ):
        super(MaxdisplayedValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
