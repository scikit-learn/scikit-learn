import _plotly_utils.basevalidators


class FeatureidkeyValidator(_plotly_utils.basevalidators.StringValidator):
    def __init__(self, plotly_name="featureidkey", parent_name="choropleth", **kwargs):
        super(FeatureidkeyValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
