import _plotly_utils.basevalidators


class TracerefminusValidator(_plotly_utils.basevalidators.IntegerValidator):
    def __init__(
        self, plotly_name="tracerefminus", parent_name="scatter3d.error_y", **kwargs
    ):
        super(TracerefminusValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
