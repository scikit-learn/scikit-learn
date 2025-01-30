import _plotly_utils.basevalidators


class BValidator(_plotly_utils.basevalidators.DataArrayValidator):
    def __init__(self, plotly_name="b", parent_name="contourcarpet", **kwargs):
        super(BValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc+clearAxisTypes"),
            implied_edits=kwargs.pop("implied_edits", {"ytype": "array"}),
            **kwargs,
        )
