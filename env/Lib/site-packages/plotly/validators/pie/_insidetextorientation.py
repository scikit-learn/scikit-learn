import _plotly_utils.basevalidators


class InsidetextorientationValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(
        self, plotly_name="insidetextorientation", parent_name="pie", **kwargs
    ):
        super(InsidetextorientationValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            values=kwargs.pop("values", ["horizontal", "radial", "tangential", "auto"]),
            **kwargs,
        )
