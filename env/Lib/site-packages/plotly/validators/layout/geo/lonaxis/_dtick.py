import _plotly_utils.basevalidators


class DtickValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="dtick", parent_name="layout.geo.lonaxis", **kwargs):
        super(DtickValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            **kwargs,
        )
