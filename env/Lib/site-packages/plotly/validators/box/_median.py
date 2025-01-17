import _plotly_utils.basevalidators


class MedianValidator(_plotly_utils.basevalidators.DataArrayValidator):
    def __init__(self, plotly_name="median", parent_name="box", **kwargs):
        super(MedianValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc+clearAxisTypes"),
            **kwargs,
        )
