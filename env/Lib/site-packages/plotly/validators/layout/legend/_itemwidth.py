import _plotly_utils.basevalidators


class ItemwidthValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="itemwidth", parent_name="layout.legend", **kwargs):
        super(ItemwidthValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "legend"),
            min=kwargs.pop("min", 30),
            **kwargs,
        )
