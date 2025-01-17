import _plotly_utils.basevalidators


class MinorgridcountValidator(_plotly_utils.basevalidators.IntegerValidator):
    def __init__(
        self, plotly_name="minorgridcount", parent_name="carpet.baxis", **kwargs
    ):
        super(MinorgridcountValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
