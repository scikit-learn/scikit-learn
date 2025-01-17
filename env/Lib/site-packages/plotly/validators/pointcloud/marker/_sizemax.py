import _plotly_utils.basevalidators


class SizemaxValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(
        self, plotly_name="sizemax", parent_name="pointcloud.marker", **kwargs
    ):
        super(SizemaxValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            min=kwargs.pop("min", 0.1),
            **kwargs,
        )
