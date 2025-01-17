import _plotly_utils.basevalidators


class YshiftValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(
        self, plotly_name="yshift", parent_name="layout.scene.annotation", **kwargs
    ):
        super(YshiftValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
