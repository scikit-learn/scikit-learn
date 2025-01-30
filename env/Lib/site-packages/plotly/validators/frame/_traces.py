import _plotly_utils.basevalidators


class TracesValidator(_plotly_utils.basevalidators.AnyValidator):
    def __init__(self, plotly_name="traces", parent_name="frame", **kwargs):
        super(TracesValidator, self).__init__(
            plotly_name=plotly_name, parent_name=parent_name, **kwargs
        )
