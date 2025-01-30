import _plotly_utils.basevalidators


class StyleValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(
        self,
        plotly_name="style",
        parent_name="scatterternary.marker.colorbar.tickfont",
        **kwargs,
    ):
        super(StyleValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "colorbars"),
            values=kwargs.pop("values", ["normal", "italic"]),
            **kwargs,
        )
