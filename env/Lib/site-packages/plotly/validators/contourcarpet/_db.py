import _plotly_utils.basevalidators


class DbValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="db", parent_name="contourcarpet", **kwargs):
        super(DbValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            implied_edits=kwargs.pop("implied_edits", {"ytype": "scaled"}),
            **kwargs,
        )
