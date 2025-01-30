import _plotly_utils.basevalidators


class DvalueValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(
        self, plotly_name="dvalue", parent_name="layout.yaxis.rangebreak", **kwargs
    ):
        super(DvalueValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
