import _plotly_utils.basevalidators


class HovercolorValidator(_plotly_utils.basevalidators.ColorValidator):
    def __init__(self, plotly_name="hovercolor", parent_name="sankey.link", **kwargs):
        super(HovercolorValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            array_ok=kwargs.pop("array_ok", True),
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
