import _plotly_utils.basevalidators


class MediansrcValidator(_plotly_utils.basevalidators.SrcValidator):
    def __init__(self, plotly_name="mediansrc", parent_name="box", **kwargs):
        super(MediansrcValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            **kwargs,
        )
