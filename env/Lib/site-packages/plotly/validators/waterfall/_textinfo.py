import _plotly_utils.basevalidators


class TextinfoValidator(_plotly_utils.basevalidators.FlaglistValidator):
    def __init__(self, plotly_name="textinfo", parent_name="waterfall", **kwargs):
        super(TextinfoValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            array_ok=kwargs.pop("array_ok", False),
            edit_type=kwargs.pop("edit_type", "plot"),
            extras=kwargs.pop("extras", ["none"]),
            flags=kwargs.pop("flags", ["label", "text", "initial", "delta", "final"]),
            **kwargs,
        )
