import _plotly_utils.basevalidators


class LegendValidator(_plotly_utils.basevalidators.SubplotidValidator):
    def __init__(self, plotly_name="legend", parent_name="scatterpolar", **kwargs):
        super(LegendValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            dflt=kwargs.pop("dflt", "legend"),
            edit_type=kwargs.pop("edit_type", "style"),
            **kwargs,
        )
