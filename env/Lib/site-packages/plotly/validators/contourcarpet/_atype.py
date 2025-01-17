import _plotly_utils.basevalidators


class AtypeValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="atype", parent_name="contourcarpet", **kwargs):
        super(AtypeValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc+clearAxisTypes"),
            values=kwargs.pop("values", ["array", "scaled"]),
            **kwargs,
        )
