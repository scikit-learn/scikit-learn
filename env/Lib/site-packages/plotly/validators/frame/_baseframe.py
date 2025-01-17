import _plotly_utils.basevalidators


class BaseframeValidator(_plotly_utils.basevalidators.StringValidator):
    def __init__(self, plotly_name="baseframe", parent_name="frame", **kwargs):
        super(BaseframeValidator, self).__init__(
            plotly_name=plotly_name, parent_name=parent_name, **kwargs
        )
