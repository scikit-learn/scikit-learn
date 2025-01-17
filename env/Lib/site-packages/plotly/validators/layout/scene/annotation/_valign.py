import _plotly_utils.basevalidators


class ValignValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(
        self, plotly_name="valign", parent_name="layout.scene.annotation", **kwargs
    ):
        super(ValignValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            values=kwargs.pop("values", ["top", "middle", "bottom"]),
            **kwargs,
        )
