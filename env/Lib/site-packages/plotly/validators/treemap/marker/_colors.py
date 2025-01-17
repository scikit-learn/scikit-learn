import _plotly_utils.basevalidators


class ColorsValidator(_plotly_utils.basevalidators.DataArrayValidator):
    def __init__(self, plotly_name="colors", parent_name="treemap.marker", **kwargs):
        super(ColorsValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
