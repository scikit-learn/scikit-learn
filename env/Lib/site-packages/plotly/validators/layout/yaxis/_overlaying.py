import _plotly_utils.basevalidators


class OverlayingValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="overlaying", parent_name="layout.yaxis", **kwargs):
        super(OverlayingValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            values=kwargs.pop(
                "values",
                [
                    "free",
                    "/^x([2-9]|[1-9][0-9]+)?( domain)?$/",
                    "/^y([2-9]|[1-9][0-9]+)?( domain)?$/",
                ],
            ),
            **kwargs,
        )
