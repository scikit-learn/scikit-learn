import _plotly_utils.basevalidators


class EndValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(
        self, plotly_name="end", parent_name="histogram2dcontour.contours", **kwargs
    ):
        super(EndValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            implied_edits=kwargs.pop("implied_edits", {"^autocontour": False}),
            **kwargs,
        )
