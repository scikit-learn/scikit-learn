import _plotly_utils.basevalidators


class ZhoverformatValidator(_plotly_utils.basevalidators.StringValidator):
    def __init__(self, plotly_name="zhoverformat", parent_name="volume", **kwargs):
        super(ZhoverformatValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
