import _plotly_utils.basevalidators


class ZsmoothValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="zsmooth", parent_name="image", **kwargs):
        super(ZsmoothValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            values=kwargs.pop("values", ["fast", False]),
            **kwargs,
        )
