import _plotly_utils.basevalidators


class SymbolValidator(_plotly_utils.basevalidators.StringValidator):
    def __init__(self, plotly_name="symbol", parent_name="scattermap.marker", **kwargs):
        super(SymbolValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            array_ok=kwargs.pop("array_ok", True),
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
