import _plotly_utils.basevalidators


class RemoveValidator(_plotly_utils.basevalidators.StringValidator):
    def __init__(self, plotly_name="remove", parent_name="layout.modebar", **kwargs):
        super(RemoveValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            array_ok=kwargs.pop("array_ok", True),
            edit_type=kwargs.pop("edit_type", "modebar"),
            **kwargs,
        )
