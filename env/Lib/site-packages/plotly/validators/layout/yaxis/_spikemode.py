import _plotly_utils.basevalidators


class SpikemodeValidator(_plotly_utils.basevalidators.FlaglistValidator):
    def __init__(self, plotly_name="spikemode", parent_name="layout.yaxis", **kwargs):
        super(SpikemodeValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            flags=kwargs.pop("flags", ["toaxis", "across", "marker"]),
            **kwargs,
        )
