import _plotly_utils.basevalidators


class XhoverformatValidator(_plotly_utils.basevalidators.StringValidator):
    def __init__(self, plotly_name="xhoverformat", parent_name="volume", **kwargs):
        super(XhoverformatValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
