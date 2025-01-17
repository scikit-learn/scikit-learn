import _plotly_utils.basevalidators


class FramewidthValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="framewidth", parent_name="layout.geo", **kwargs):
        super(FramewidthValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
