import _plotly_utils.basevalidators


class AyValidator(_plotly_utils.basevalidators.AnyValidator):
    def __init__(self, plotly_name="ay", parent_name="layout.annotation", **kwargs):
        super(AyValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc+arraydraw"),
            **kwargs,
        )
