import _plotly_utils.basevalidators


class AmbientValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="ambient", parent_name="volume.lighting", **kwargs):
        super(AmbientValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            max=kwargs.pop("max", 1),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
