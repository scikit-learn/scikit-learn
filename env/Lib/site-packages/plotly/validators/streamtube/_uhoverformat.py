import _plotly_utils.basevalidators


class UhoverformatValidator(_plotly_utils.basevalidators.StringValidator):
    def __init__(self, plotly_name="uhoverformat", parent_name="streamtube", **kwargs):
        super(UhoverformatValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            **kwargs,
        )
