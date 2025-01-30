import _plotly_utils.basevalidators


class XaxesValidator(_plotly_utils.basevalidators.InfoArrayValidator):
    def __init__(self, plotly_name="xaxes", parent_name="splom", **kwargs):
        super(XaxesValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            free_length=kwargs.pop("free_length", True),
            items=kwargs.pop(
                "items",
                {
                    "editType": "plot",
                    "regex": "/^x([2-9]|[1-9][0-9]+)?( domain)?$/",
                    "valType": "subplotid",
                },
            ),
            **kwargs,
        )
