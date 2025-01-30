import _plotly_utils.basevalidators


class IntensitymodeValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="intensitymode", parent_name="mesh3d", **kwargs):
        super(IntensitymodeValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            values=kwargs.pop("values", ["vertex", "cell"]),
            **kwargs,
        )
