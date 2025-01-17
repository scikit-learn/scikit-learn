import _plotly_utils.basevalidators


class LabelValidator(_plotly_utils.basevalidators.StringValidator):
    def __init__(self, plotly_name="label", parent_name="layout.slider.step", **kwargs):
        super(LabelValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "arraydraw"),
            **kwargs,
        )
