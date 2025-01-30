import _plotly_utils.basevalidators


class SizerefValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="sizeref", parent_name="streamtube", **kwargs):
        super(SizerefValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
