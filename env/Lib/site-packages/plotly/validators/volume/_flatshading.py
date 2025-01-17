import _plotly_utils.basevalidators


class FlatshadingValidator(_plotly_utils.basevalidators.BooleanValidator):
    def __init__(self, plotly_name="flatshading", parent_name="volume", **kwargs):
        super(FlatshadingValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
