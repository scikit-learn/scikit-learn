import plotly.validators


class LayoutValidator(plotly.validators.LayoutValidator):
    def __init__(self, plotly_name="layout", parent_name="frame", **kwargs):
        super(LayoutValidator, self).__init__(
            plotly_name=plotly_name, parent_name=parent_name, **kwargs
        )
