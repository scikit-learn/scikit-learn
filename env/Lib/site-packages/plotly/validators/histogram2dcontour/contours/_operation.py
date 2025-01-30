import _plotly_utils.basevalidators


class OperationValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(
        self,
        plotly_name="operation",
        parent_name="histogram2dcontour.contours",
        **kwargs,
    ):
        super(OperationValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            values=kwargs.pop(
                "values",
                [
                    "=",
                    "<",
                    ">=",
                    ">",
                    "<=",
                    "[]",
                    "()",
                    "[)",
                    "(]",
                    "][",
                    ")(",
                    "](",
                    ")[",
                ],
            ),
            **kwargs,
        )
