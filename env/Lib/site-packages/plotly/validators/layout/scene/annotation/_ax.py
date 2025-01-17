import _plotly_utils.basevalidators


class AxValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(
        self, plotly_name="ax", parent_name="layout.scene.annotation", **kwargs
    ):
        super(AxValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
