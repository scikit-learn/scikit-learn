import _plotly_utils.basevalidators


class Theta0Validator(_plotly_utils.basevalidators.AnyValidator):
    def __init__(self, plotly_name="theta0", parent_name="barpolar", **kwargs):
        super(Theta0Validator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc+clearAxisTypes"),
            **kwargs,
        )
