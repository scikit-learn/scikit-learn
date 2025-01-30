import _plotly_utils.basevalidators


class ShowlinesValidator(_plotly_utils.basevalidators.BooleanValidator):
    def __init__(
        self,
        plotly_name="showlines",
        parent_name="histogram2dcontour.contours",
        **kwargs,
    ):
        super(ShowlinesValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            **kwargs,
        )
