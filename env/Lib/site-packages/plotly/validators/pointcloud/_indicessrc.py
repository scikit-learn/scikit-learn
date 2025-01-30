import _plotly_utils.basevalidators


class IndicessrcValidator(_plotly_utils.basevalidators.SrcValidator):
    def __init__(self, plotly_name="indicessrc", parent_name="pointcloud", **kwargs):
        super(IndicessrcValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            **kwargs,
        )
