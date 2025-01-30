import _plotly_utils.basevalidators


class ShadowValidator(_plotly_utils.basevalidators.StringValidator):
    def __init__(self, plotly_name="shadow", parent_name="table.header.font", **kwargs):
        super(ShadowValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            array_ok=kwargs.pop("array_ok", True),
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
