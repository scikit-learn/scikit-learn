import _plotly_utils.basevalidators


class CornerradiusValidator(_plotly_utils.basevalidators.AnyValidator):
    def __init__(self, plotly_name="cornerradius", parent_name="bar.marker", **kwargs):
        super(CornerradiusValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
