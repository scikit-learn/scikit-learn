import _plotly_utils.basevalidators


class ArraydtickValidator(_plotly_utils.basevalidators.IntegerValidator):
    def __init__(self, plotly_name="arraydtick", parent_name="carpet.baxis", **kwargs):
        super(ArraydtickValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            min=kwargs.pop("min", 1),
            **kwargs,
        )
