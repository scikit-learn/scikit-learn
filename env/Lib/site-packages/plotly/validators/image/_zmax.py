import _plotly_utils.basevalidators


class ZmaxValidator(_plotly_utils.basevalidators.InfoArrayValidator):
    def __init__(self, plotly_name="zmax", parent_name="image", **kwargs):
        super(ZmaxValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            items=kwargs.pop(
                "items",
                [
                    {"editType": "calc", "valType": "number"},
                    {"editType": "calc", "valType": "number"},
                    {"editType": "calc", "valType": "number"},
                    {"editType": "calc", "valType": "number"},
                ],
            ),
            **kwargs,
        )
