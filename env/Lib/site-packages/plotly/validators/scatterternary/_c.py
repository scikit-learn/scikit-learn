import _plotly_utils.basevalidators


class CValidator(_plotly_utils.basevalidators.DataArrayValidator):
    def __init__(self, plotly_name="c", parent_name="scatterternary", **kwargs):
        super(CValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
