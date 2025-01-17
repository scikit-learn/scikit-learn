import _plotly_utils.basevalidators


class ResolutionValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="resolution", parent_name="layout.geo", **kwargs):
        super(ResolutionValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            coerce_number=kwargs.pop("coerce_number", True),
            edit_type=kwargs.pop("edit_type", "plot"),
            values=kwargs.pop("values", [110, 50]),
            **kwargs,
        )
