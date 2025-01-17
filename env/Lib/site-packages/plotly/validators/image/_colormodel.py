import _plotly_utils.basevalidators


class ColormodelValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="colormodel", parent_name="image", **kwargs):
        super(ColormodelValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            values=kwargs.pop("values", ["rgb", "rgba", "rgba256", "hsl", "hsla"]),
            **kwargs,
        )
