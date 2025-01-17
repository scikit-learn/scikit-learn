import _plotly_utils.basevalidators


class SizerefValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="sizeref", parent_name="splom.marker", **kwargs):
        super(SizerefValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
