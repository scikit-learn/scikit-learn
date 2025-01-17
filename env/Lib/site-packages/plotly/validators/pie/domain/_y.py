import _plotly_utils.basevalidators


class YValidator(_plotly_utils.basevalidators.InfoArrayValidator):
    def __init__(self, plotly_name="y", parent_name="pie.domain", **kwargs):
        super(YValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            items=kwargs.pop(
                "items",
                [
                    {"editType": "calc", "max": 1, "min": 0, "valType": "number"},
                    {"editType": "calc", "max": 1, "min": 0, "valType": "number"},
                ],
            ),
            **kwargs,
        )
