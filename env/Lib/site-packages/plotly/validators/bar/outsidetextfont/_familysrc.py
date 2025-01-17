import _plotly_utils.basevalidators


class FamilysrcValidator(_plotly_utils.basevalidators.SrcValidator):
    def __init__(
        self, plotly_name="familysrc", parent_name="bar.outsidetextfont", **kwargs
    ):
        super(FamilysrcValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            **kwargs,
        )
