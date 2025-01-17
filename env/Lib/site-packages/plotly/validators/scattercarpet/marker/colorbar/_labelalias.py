import _plotly_utils.basevalidators


class LabelaliasValidator(_plotly_utils.basevalidators.AnyValidator):
    def __init__(
        self,
        plotly_name="labelalias",
        parent_name="scattercarpet.marker.colorbar",
        **kwargs,
    ):
        super(LabelaliasValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "colorbars"),
            **kwargs,
        )
