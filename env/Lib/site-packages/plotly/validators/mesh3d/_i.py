import _plotly_utils.basevalidators


class IValidator(_plotly_utils.basevalidators.DataArrayValidator):
    def __init__(self, plotly_name="i", parent_name="mesh3d", **kwargs):
        super(IValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
