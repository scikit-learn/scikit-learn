import _plotly_utils.basevalidators


class MinorticklenValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(
        self, plotly_name="minorticklen", parent_name="layout.slider", **kwargs
    ):
        super(MinorticklenValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "arraydraw"),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
