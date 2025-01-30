import _plotly_utils.basevalidators


class ShowarrowValidator(_plotly_utils.basevalidators.BooleanValidator):
    def __init__(
        self, plotly_name="showarrow", parent_name="layout.scene.annotation", **kwargs
    ):
        super(ShowarrowValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
