import _plotly_utils.basevalidators


class ThetaValidator(_plotly_utils.basevalidators.DataArrayValidator):
    def __init__(self, plotly_name="theta", parent_name="scatterpolar", **kwargs):
        super(ThetaValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc+clearAxisTypes"),
            **kwargs,
        )
