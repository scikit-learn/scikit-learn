import _plotly_utils.basevalidators


class EdgeshapeValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="edgeshape", parent_name="icicle.pathbar", **kwargs):
        super(EdgeshapeValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            values=kwargs.pop("values", [">", "<", "|", "/", "\\"]),
            **kwargs,
        )
