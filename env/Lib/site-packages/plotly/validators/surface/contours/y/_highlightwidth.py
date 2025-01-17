import _plotly_utils.basevalidators


class HighlightwidthValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(
        self, plotly_name="highlightwidth", parent_name="surface.contours.y", **kwargs
    ):
        super(HighlightwidthValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            max=kwargs.pop("max", 16),
            min=kwargs.pop("min", 1),
            **kwargs,
        )
