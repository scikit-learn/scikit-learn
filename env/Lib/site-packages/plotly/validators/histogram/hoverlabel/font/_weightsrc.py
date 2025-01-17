import _plotly_utils.basevalidators


class WeightsrcValidator(_plotly_utils.basevalidators.SrcValidator):
    def __init__(
        self, plotly_name="weightsrc", parent_name="histogram.hoverlabel.font", **kwargs
    ):
        super(WeightsrcValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            **kwargs,
        )
