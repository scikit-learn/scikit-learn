import _plotly_utils.basevalidators


class B0Validator(_plotly_utils.basevalidators.AnyValidator):
    def __init__(self, plotly_name="b0", parent_name="contourcarpet", **kwargs):
        super(B0Validator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc+clearAxisTypes"),
            implied_edits=kwargs.pop("implied_edits", {"ytype": "scaled"}),
            **kwargs,
        )
