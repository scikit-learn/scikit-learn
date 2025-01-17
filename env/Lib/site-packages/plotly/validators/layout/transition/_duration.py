import _plotly_utils.basevalidators


class DurationValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(
        self, plotly_name="duration", parent_name="layout.transition", **kwargs
    ):
        super(DurationValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
