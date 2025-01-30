import _plotly_utils.basevalidators


class AyrefValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="ayref", parent_name="layout.annotation", **kwargs):
        super(AyrefValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            values=kwargs.pop(
                "values", ["pixel", "/^y([2-9]|[1-9][0-9]+)?( domain)?$/"]
            ),
            **kwargs,
        )
