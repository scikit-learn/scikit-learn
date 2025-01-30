import _plotly_utils.basevalidators


class DxValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="dx", parent_name="bar", **kwargs):
        super(DxValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            anim=kwargs.pop("anim", True),
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
