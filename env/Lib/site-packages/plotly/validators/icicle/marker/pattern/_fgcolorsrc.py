import _plotly_utils.basevalidators


class FgcolorsrcValidator(_plotly_utils.basevalidators.SrcValidator):
    def __init__(
        self, plotly_name="fgcolorsrc", parent_name="icicle.marker.pattern", **kwargs
    ):
        super(FgcolorsrcValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            **kwargs,
        )
