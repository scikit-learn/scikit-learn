import _plotly_utils.basevalidators


class ConstraintowardValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(
        self, plotly_name="constraintoward", parent_name="layout.yaxis", **kwargs
    ):
        super(ConstraintowardValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            values=kwargs.pop(
                "values", ["left", "center", "right", "top", "middle", "bottom"]
            ),
            **kwargs,
        )
