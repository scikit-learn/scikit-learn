import _plotly_utils.basevalidators


class ClickmodeValidator(_plotly_utils.basevalidators.FlaglistValidator):
    def __init__(self, plotly_name="clickmode", parent_name="layout", **kwargs):
        super(ClickmodeValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            extras=kwargs.pop("extras", ["none"]),
            flags=kwargs.pop("flags", ["event", "select"]),
            **kwargs,
        )
