import _plotly_utils.basevalidators


class YhoverformatValidator(_plotly_utils.basevalidators.StringValidator):
    def __init__(
        self, plotly_name="yhoverformat", parent_name="histogram2dcontour", **kwargs
    ):
        super(YhoverformatValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            **kwargs,
        )
