import _plotly_utils.basevalidators


class ArraysrcValidator(_plotly_utils.basevalidators.SrcValidator):
    def __init__(self, plotly_name="arraysrc", parent_name="bar.error_y", **kwargs):
        super(ArraysrcValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            **kwargs,
        )
