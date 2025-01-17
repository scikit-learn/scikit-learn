import _plotly_utils.basevalidators


class PaddingValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(
        self, plotly_name="padding", parent_name="layout.newshape.label", **kwargs
    ):
        super(PaddingValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
