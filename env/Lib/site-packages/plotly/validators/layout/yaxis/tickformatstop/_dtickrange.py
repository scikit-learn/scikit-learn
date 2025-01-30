import _plotly_utils.basevalidators


class DtickrangeValidator(_plotly_utils.basevalidators.InfoArrayValidator):
    def __init__(
        self,
        plotly_name="dtickrange",
        parent_name="layout.yaxis.tickformatstop",
        **kwargs,
    ):
        super(DtickrangeValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "ticks"),
            items=kwargs.pop(
                "items",
                [
                    {"editType": "ticks", "valType": "any"},
                    {"editType": "ticks", "valType": "any"},
                ],
            ),
            **kwargs,
        )
