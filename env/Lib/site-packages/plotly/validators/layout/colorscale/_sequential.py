import _plotly_utils.basevalidators


class SequentialValidator(_plotly_utils.basevalidators.ColorscaleValidator):
    def __init__(
        self, plotly_name="sequential", parent_name="layout.colorscale", **kwargs
    ):
        super(SequentialValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
