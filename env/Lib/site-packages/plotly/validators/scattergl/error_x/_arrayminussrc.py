import _plotly_utils.basevalidators


class ArrayminussrcValidator(_plotly_utils.basevalidators.SrcValidator):
    def __init__(
        self, plotly_name="arrayminussrc", parent_name="scattergl.error_x", **kwargs
    ):
        super(ArrayminussrcValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            **kwargs,
        )
