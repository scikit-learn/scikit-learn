import _plotly_utils.basevalidators


class AutomarginValidator(_plotly_utils.basevalidators.FlaglistValidator):
    def __init__(self, plotly_name="automargin", parent_name="layout.yaxis", **kwargs):
        super(AutomarginValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "ticks"),
            extras=kwargs.pop("extras", [True, False]),
            flags=kwargs.pop(
                "flags", ["height", "width", "left", "right", "top", "bottom"]
            ),
            **kwargs,
        )
