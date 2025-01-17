import _plotly_utils.basevalidators


class ConstrainValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="constrain", parent_name="layout.xaxis", **kwargs):
        super(ConstrainValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            values=kwargs.pop("values", ["range", "domain"]),
            **kwargs,
        )
