import _plotly_utils.basevalidators


class CmidValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="cmid", parent_name="scatter3d.marker", **kwargs):
        super(CmidValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            implied_edits=kwargs.pop("implied_edits", {}),
            **kwargs,
        )
