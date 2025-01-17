import _plotly_utils.basevalidators


class ScattermodeValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="scattermode", parent_name="layout", **kwargs):
        super(ScattermodeValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            values=kwargs.pop("values", ["group", "overlay"]),
            **kwargs,
        )
