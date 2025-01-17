import _plotly_utils.basevalidators


class HoveronValidator(_plotly_utils.basevalidators.FlaglistValidator):
    def __init__(self, plotly_name="hoveron", parent_name="scatterternary", **kwargs):
        super(HoveronValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "style"),
            flags=kwargs.pop("flags", ["points", "fills"]),
            **kwargs,
        )
