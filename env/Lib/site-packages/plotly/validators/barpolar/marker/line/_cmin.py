import _plotly_utils.basevalidators


class CminValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(
        self, plotly_name="cmin", parent_name="barpolar.marker.line", **kwargs
    ):
        super(CminValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            implied_edits=kwargs.pop("implied_edits", {"cauto": False}),
            **kwargs,
        )
