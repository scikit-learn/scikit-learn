import _plotly_utils.basevalidators


class LabelformatValidator(_plotly_utils.basevalidators.StringValidator):
    def __init__(
        self,
        plotly_name="labelformat",
        parent_name="histogram2dcontour.contours",
        **kwargs,
    ):
        super(LabelformatValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            **kwargs,
        )
