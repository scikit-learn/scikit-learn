import _plotly_utils.basevalidators


class TraceorderValidator(_plotly_utils.basevalidators.FlaglistValidator):
    def __init__(self, plotly_name="traceorder", parent_name="layout.legend", **kwargs):
        super(TraceorderValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "legend"),
            extras=kwargs.pop("extras", ["normal"]),
            flags=kwargs.pop("flags", ["reversed", "grouped"]),
            **kwargs,
        )
