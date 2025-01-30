import _plotly_utils.basevalidators


class DlabelValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="dlabel", parent_name="funnelarea", **kwargs):
        super(DlabelValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
