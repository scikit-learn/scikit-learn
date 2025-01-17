import _plotly_utils.basevalidators


class ShowdividersValidator(_plotly_utils.basevalidators.BooleanValidator):
    def __init__(
        self, plotly_name="showdividers", parent_name="layout.yaxis", **kwargs
    ):
        super(ShowdividersValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "ticks"),
            **kwargs,
        )
