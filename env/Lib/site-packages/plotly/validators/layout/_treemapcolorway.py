import _plotly_utils.basevalidators


class TreemapcolorwayValidator(_plotly_utils.basevalidators.ColorlistValidator):
    def __init__(self, plotly_name="treemapcolorway", parent_name="layout", **kwargs):
        super(TreemapcolorwayValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
