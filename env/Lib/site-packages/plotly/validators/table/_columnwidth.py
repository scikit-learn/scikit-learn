import _plotly_utils.basevalidators


class ColumnwidthValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="columnwidth", parent_name="table", **kwargs):
        super(ColumnwidthValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            array_ok=kwargs.pop("array_ok", True),
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
