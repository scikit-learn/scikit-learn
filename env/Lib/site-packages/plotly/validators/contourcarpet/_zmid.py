import _plotly_utils.basevalidators


class ZmidValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="zmid", parent_name="contourcarpet", **kwargs):
        super(ZmidValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            implied_edits=kwargs.pop("implied_edits", {}),
            **kwargs,
        )
