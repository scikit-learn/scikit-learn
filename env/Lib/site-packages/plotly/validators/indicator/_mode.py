import _plotly_utils.basevalidators


class ModeValidator(_plotly_utils.basevalidators.FlaglistValidator):
    def __init__(self, plotly_name="mode", parent_name="indicator", **kwargs):
        super(ModeValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            flags=kwargs.pop("flags", ["number", "delta", "gauge"]),
            **kwargs,
        )
