import _plotly_utils.basevalidators


class AngleValidator(_plotly_utils.basevalidators.AngleValidator):
    def __init__(
        self, plotly_name="angle", parent_name="layout.polar.radialaxis", **kwargs
    ):
        super(AngleValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            **kwargs,
        )
