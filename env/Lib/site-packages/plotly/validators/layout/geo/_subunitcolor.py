import _plotly_utils.basevalidators


class SubunitcolorValidator(_plotly_utils.basevalidators.ColorValidator):
    def __init__(self, plotly_name="subunitcolor", parent_name="layout.geo", **kwargs):
        super(SubunitcolorValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            **kwargs,
        )
