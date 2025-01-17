import _plotly_utils.basevalidators


class WidthsrcValidator(_plotly_utils.basevalidators.SrcValidator):
    def __init__(
        self, plotly_name="widthsrc", parent_name="funnel.marker.line", **kwargs
    ):
        super(WidthsrcValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            **kwargs,
        )
