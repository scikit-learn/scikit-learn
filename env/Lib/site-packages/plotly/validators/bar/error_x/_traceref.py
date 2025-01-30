import _plotly_utils.basevalidators


class TracerefValidator(_plotly_utils.basevalidators.IntegerValidator):
    def __init__(self, plotly_name="traceref", parent_name="bar.error_x", **kwargs):
        super(TracerefValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "style"),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
