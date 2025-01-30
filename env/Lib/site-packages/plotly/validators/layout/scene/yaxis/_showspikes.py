import _plotly_utils.basevalidators


class ShowspikesValidator(_plotly_utils.basevalidators.BooleanValidator):
    def __init__(
        self, plotly_name="showspikes", parent_name="layout.scene.yaxis", **kwargs
    ):
        super(ShowspikesValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            **kwargs,
        )
