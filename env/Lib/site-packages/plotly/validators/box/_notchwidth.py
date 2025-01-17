import _plotly_utils.basevalidators


class NotchwidthValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="notchwidth", parent_name="box", **kwargs):
        super(NotchwidthValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            max=kwargs.pop("max", 0.5),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
