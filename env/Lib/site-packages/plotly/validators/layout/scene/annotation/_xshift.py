import _plotly_utils.basevalidators


class XshiftValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(
        self, plotly_name="xshift", parent_name="layout.scene.annotation", **kwargs
    ):
        super(XshiftValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
