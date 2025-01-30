import _plotly_utils.basevalidators


class IsomaxValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="isomax", parent_name="isosurface", **kwargs):
        super(IsomaxValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
