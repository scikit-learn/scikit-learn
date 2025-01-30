import _plotly_utils.basevalidators


class HovermodeValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="hovermode", parent_name="layout", **kwargs):
        super(HovermodeValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "modebar"),
            values=kwargs.pop(
                "values", ["x", "y", "closest", False, "x unified", "y unified"]
            ),
            **kwargs,
        )
