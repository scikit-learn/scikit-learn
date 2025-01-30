import _plotly_utils.basevalidators


class DisplayindexValidator(_plotly_utils.basevalidators.IntegerValidator):
    def __init__(
        self, plotly_name="displayindex", parent_name="parcats.dimension", **kwargs
    ):
        super(DisplayindexValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
