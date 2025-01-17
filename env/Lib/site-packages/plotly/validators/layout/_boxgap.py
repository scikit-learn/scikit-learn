import _plotly_utils.basevalidators


class BoxgapValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="boxgap", parent_name="layout", **kwargs):
        super(BoxgapValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            max=kwargs.pop("max", 1),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
