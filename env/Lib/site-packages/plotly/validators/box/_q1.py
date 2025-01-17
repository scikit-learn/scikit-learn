import _plotly_utils.basevalidators


class Q1Validator(_plotly_utils.basevalidators.DataArrayValidator):
    def __init__(self, plotly_name="q1", parent_name="box", **kwargs):
        super(Q1Validator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc+clearAxisTypes"),
            **kwargs,
        )
