import _plotly_utils.basevalidators


class NamelengthValidator(_plotly_utils.basevalidators.IntegerValidator):
    def __init__(
        self, plotly_name="namelength", parent_name="ohlc.hoverlabel", **kwargs
    ):
        super(NamelengthValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            array_ok=kwargs.pop("array_ok", True),
            edit_type=kwargs.pop("edit_type", "none"),
            min=kwargs.pop("min", -1),
            **kwargs,
        )
