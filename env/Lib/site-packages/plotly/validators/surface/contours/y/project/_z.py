import _plotly_utils.basevalidators


class ZValidator(_plotly_utils.basevalidators.BooleanValidator):
    def __init__(
        self, plotly_name="z", parent_name="surface.contours.y.project", **kwargs
    ):
        super(ZValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
