import _plotly_utils.basevalidators


class AxrefValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="axref", parent_name="layout.annotation", **kwargs):
        super(AxrefValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            values=kwargs.pop(
                "values", ["pixel", "/^x([2-9]|[1-9][0-9]+)?( domain)?$/"]
            ),
            **kwargs,
        )
