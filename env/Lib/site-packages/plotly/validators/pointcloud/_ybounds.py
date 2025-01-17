import _plotly_utils.basevalidators


class YboundsValidator(_plotly_utils.basevalidators.DataArrayValidator):
    def __init__(self, plotly_name="ybounds", parent_name="pointcloud", **kwargs):
        super(YboundsValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
