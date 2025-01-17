import _plotly_utils.basevalidators


class ClicktoshowValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(
        self, plotly_name="clicktoshow", parent_name="layout.annotation", **kwargs
    ):
        super(ClicktoshowValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "arraydraw"),
            values=kwargs.pop("values", [False, "onoff", "onout"]),
            **kwargs,
        )
