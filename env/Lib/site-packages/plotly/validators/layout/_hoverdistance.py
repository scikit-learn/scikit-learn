import _plotly_utils.basevalidators


class HoverdistanceValidator(_plotly_utils.basevalidators.IntegerValidator):
    def __init__(self, plotly_name="hoverdistance", parent_name="layout", **kwargs):
        super(HoverdistanceValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            min=kwargs.pop("min", -1),
            **kwargs,
        )
