import _plotly_utils.basevalidators


class TextangleValidator(_plotly_utils.basevalidators.AngleValidator):
    def __init__(
        self, plotly_name="textangle", parent_name="layout.scene.annotation", **kwargs
    ):
        super(TextangleValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
