import _plotly_utils.basevalidators


class ThicknessmodeValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(
        self,
        plotly_name="thicknessmode",
        parent_name="choroplethmapbox.colorbar",
        **kwargs,
    ):
        super(ThicknessmodeValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "colorbars"),
            values=kwargs.pop("values", ["fraction", "pixels"]),
            **kwargs,
        )
