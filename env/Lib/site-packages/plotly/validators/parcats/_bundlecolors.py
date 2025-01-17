import _plotly_utils.basevalidators


class BundlecolorsValidator(_plotly_utils.basevalidators.BooleanValidator):
    def __init__(self, plotly_name="bundlecolors", parent_name="parcats", **kwargs):
        super(BundlecolorsValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            **kwargs,
        )
