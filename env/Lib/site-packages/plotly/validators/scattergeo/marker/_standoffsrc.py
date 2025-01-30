import _plotly_utils.basevalidators


class StandoffsrcValidator(_plotly_utils.basevalidators.SrcValidator):
    def __init__(
        self, plotly_name="standoffsrc", parent_name="scattergeo.marker", **kwargs
    ):
        super(StandoffsrcValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            **kwargs,
        )
