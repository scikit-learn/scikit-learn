import _plotly_utils.basevalidators


class MinreducedwidthValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="minreducedwidth", parent_name="layout", **kwargs):
        super(MinreducedwidthValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            min=kwargs.pop("min", 2),
            **kwargs,
        )
