import _plotly_utils.basevalidators


class ColumnValidator(_plotly_utils.basevalidators.IntegerValidator):
    def __init__(self, plotly_name="column", parent_name="parcoords.domain", **kwargs):
        super(ColumnValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
