import _plotly_utils.basevalidators


class AnglesrcValidator(_plotly_utils.basevalidators.SrcValidator):
    def __init__(
        self, plotly_name="anglesrc", parent_name="scattermapbox.marker", **kwargs
    ):
        super(AnglesrcValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            **kwargs,
        )
