import _plotly_utils.basevalidators


class ZValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(
        self, plotly_name="z", parent_name="isosurface.lightposition", **kwargs
    ):
        super(ZValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            max=kwargs.pop("max", 100000),
            min=kwargs.pop("min", -100000),
            **kwargs,
        )
