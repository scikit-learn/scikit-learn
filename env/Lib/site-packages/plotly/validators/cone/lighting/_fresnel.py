import _plotly_utils.basevalidators


class FresnelValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="fresnel", parent_name="cone.lighting", **kwargs):
        super(FresnelValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            max=kwargs.pop("max", 5),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
