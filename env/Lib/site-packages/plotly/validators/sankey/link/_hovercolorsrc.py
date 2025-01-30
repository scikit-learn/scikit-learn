import _plotly_utils.basevalidators


class HovercolorsrcValidator(_plotly_utils.basevalidators.SrcValidator):
    def __init__(
        self, plotly_name="hovercolorsrc", parent_name="sankey.link", **kwargs
    ):
        super(HovercolorsrcValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            **kwargs,
        )
