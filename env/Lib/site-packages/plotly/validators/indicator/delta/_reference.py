import _plotly_utils.basevalidators


class ReferenceValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(
        self, plotly_name="reference", parent_name="indicator.delta", **kwargs
    ):
        super(ReferenceValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
