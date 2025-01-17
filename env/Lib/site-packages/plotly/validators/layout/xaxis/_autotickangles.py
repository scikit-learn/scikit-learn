import _plotly_utils.basevalidators


class AutotickanglesValidator(_plotly_utils.basevalidators.InfoArrayValidator):
    def __init__(
        self, plotly_name="autotickangles", parent_name="layout.xaxis", **kwargs
    ):
        super(AutotickanglesValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "ticks"),
            free_length=kwargs.pop("free_length", True),
            items=kwargs.pop("items", {"valType": "angle"}),
            **kwargs,
        )
