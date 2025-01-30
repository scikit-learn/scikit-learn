import _plotly_utils.basevalidators


class ThetaunitValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="thetaunit", parent_name="scatterpolar", **kwargs):
        super(ThetaunitValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc+clearAxisTypes"),
            values=kwargs.pop("values", ["radians", "degrees", "gradians"]),
            **kwargs,
        )
