import _plotly_utils.basevalidators


class SpikesnapValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="spikesnap", parent_name="layout.xaxis", **kwargs):
        super(SpikesnapValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            values=kwargs.pop("values", ["data", "cursor", "hovered data"]),
            **kwargs,
        )
