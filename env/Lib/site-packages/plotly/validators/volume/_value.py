import _plotly_utils.basevalidators


class ValueValidator(_plotly_utils.basevalidators.DataArrayValidator):
    def __init__(self, plotly_name="value", parent_name="volume", **kwargs):
        super(ValueValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc+clearAxisTypes"),
            **kwargs,
        )
