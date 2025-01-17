import _plotly_utils.basevalidators


class LocationssrcValidator(_plotly_utils.basevalidators.SrcValidator):
    def __init__(
        self, plotly_name="locationssrc", parent_name="volume.slices.y", **kwargs
    ):
        super(LocationssrcValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            **kwargs,
        )
