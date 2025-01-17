import _plotly_utils.basevalidators


class StartarrowheadValidator(_plotly_utils.basevalidators.IntegerValidator):
    def __init__(
        self, plotly_name="startarrowhead", parent_name="layout.annotation", **kwargs
    ):
        super(StartarrowheadValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "arraydraw"),
            max=kwargs.pop("max", 8),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
