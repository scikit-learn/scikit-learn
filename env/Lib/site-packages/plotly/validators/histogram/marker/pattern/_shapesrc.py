import _plotly_utils.basevalidators


class ShapesrcValidator(_plotly_utils.basevalidators.SrcValidator):
    def __init__(
        self, plotly_name="shapesrc", parent_name="histogram.marker.pattern", **kwargs
    ):
        super(ShapesrcValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            **kwargs,
        )
