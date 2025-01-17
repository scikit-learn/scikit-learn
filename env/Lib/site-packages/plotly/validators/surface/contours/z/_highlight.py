import _plotly_utils.basevalidators


class HighlightValidator(_plotly_utils.basevalidators.BooleanValidator):
    def __init__(
        self, plotly_name="highlight", parent_name="surface.contours.z", **kwargs
    ):
        super(HighlightValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
