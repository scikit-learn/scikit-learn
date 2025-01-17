import _plotly_utils.basevalidators


class FramecolorValidator(_plotly_utils.basevalidators.ColorValidator):
    def __init__(self, plotly_name="framecolor", parent_name="layout.geo", **kwargs):
        super(FramecolorValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            **kwargs,
        )
