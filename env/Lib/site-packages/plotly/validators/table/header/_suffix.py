import _plotly_utils.basevalidators


class SuffixValidator(_plotly_utils.basevalidators.StringValidator):
    def __init__(self, plotly_name="suffix", parent_name="table.header", **kwargs):
        super(SuffixValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            array_ok=kwargs.pop("array_ok", True),
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
