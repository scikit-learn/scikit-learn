import _plotly_utils.basevalidators


class LinepositionsrcValidator(_plotly_utils.basevalidators.SrcValidator):
    def __init__(
        self, plotly_name="linepositionsrc", parent_name="scatter.textfont", **kwargs
    ):
        super(LinepositionsrcValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            **kwargs,
        )
