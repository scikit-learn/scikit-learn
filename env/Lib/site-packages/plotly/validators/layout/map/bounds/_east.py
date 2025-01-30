import _plotly_utils.basevalidators


class EastValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="east", parent_name="layout.map.bounds", **kwargs):
        super(EastValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            **kwargs,
        )
