import _plotly_utils.basevalidators


class YtypeValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="ytype", parent_name="contour", **kwargs):
        super(YtypeValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc+clearAxisTypes"),
            values=kwargs.pop("values", ["array", "scaled"]),
            **kwargs,
        )
