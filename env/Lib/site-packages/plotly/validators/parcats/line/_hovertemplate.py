import _plotly_utils.basevalidators


class HovertemplateValidator(_plotly_utils.basevalidators.StringValidator):
    def __init__(
        self, plotly_name="hovertemplate", parent_name="parcats.line", **kwargs
    ):
        super(HovertemplateValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            **kwargs,
        )
