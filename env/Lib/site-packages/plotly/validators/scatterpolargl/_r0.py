import _plotly_utils.basevalidators


class R0Validator(_plotly_utils.basevalidators.AnyValidator):
    def __init__(self, plotly_name="r0", parent_name="scatterpolargl", **kwargs):
        super(R0Validator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc+clearAxisTypes"),
            **kwargs,
        )
