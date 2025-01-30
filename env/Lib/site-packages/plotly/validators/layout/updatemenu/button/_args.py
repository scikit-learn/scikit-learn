import _plotly_utils.basevalidators


class ArgsValidator(_plotly_utils.basevalidators.InfoArrayValidator):
    def __init__(
        self, plotly_name="args", parent_name="layout.updatemenu.button", **kwargs
    ):
        super(ArgsValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "arraydraw"),
            free_length=kwargs.pop("free_length", True),
            items=kwargs.pop(
                "items",
                [
                    {"editType": "arraydraw", "valType": "any"},
                    {"editType": "arraydraw", "valType": "any"},
                    {"editType": "arraydraw", "valType": "any"},
                ],
            ),
            **kwargs,
        )
