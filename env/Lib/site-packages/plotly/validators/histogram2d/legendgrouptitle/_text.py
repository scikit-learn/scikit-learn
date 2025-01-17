import _plotly_utils.basevalidators


class TextValidator(_plotly_utils.basevalidators.StringValidator):
    def __init__(
        self, plotly_name="text", parent_name="histogram2d.legendgrouptitle", **kwargs
    ):
        super(TextValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "style"),
            **kwargs,
        )
