import _plotly_utils.basevalidators


class ValuesValidator(_plotly_utils.basevalidators.InfoArrayValidator):
    def __init__(
        self, plotly_name="values", parent_name="layout.yaxis.rangebreak", **kwargs
    ):
        super(ValuesValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            free_length=kwargs.pop("free_length", True),
            items=kwargs.pop("items", {"editType": "calc", "valType": "any"}),
            **kwargs,
        )
