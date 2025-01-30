import _plotly_utils.basevalidators


class OffsetValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(
        self, plotly_name="offset", parent_name="carpet.baxis.title", **kwargs
    ):
        super(OffsetValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
