import _plotly_utils.basevalidators


class ArrowlenValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="arrowlen", parent_name="sankey.link", **kwargs):
        super(ArrowlenValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
