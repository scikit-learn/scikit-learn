import _plotly_utils.basevalidators


class RValidator(_plotly_utils.basevalidators.DataArrayValidator):
    def __init__(self, plotly_name="r", parent_name="scatterpolargl", **kwargs):
        super(RValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc+clearAxisTypes"),
            **kwargs,
        )
