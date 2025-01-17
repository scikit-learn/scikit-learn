import _plotly_utils.basevalidators


class EndlinecolorValidator(_plotly_utils.basevalidators.ColorValidator):
    def __init__(
        self, plotly_name="endlinecolor", parent_name="carpet.baxis", **kwargs
    ):
        super(EndlinecolorValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
