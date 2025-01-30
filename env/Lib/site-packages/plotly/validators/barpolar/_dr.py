import _plotly_utils.basevalidators


class DrValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="dr", parent_name="barpolar", **kwargs):
        super(DrValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
