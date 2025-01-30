import _plotly_utils.basevalidators


class SelectionrevisionValidator(_plotly_utils.basevalidators.AnyValidator):
    def __init__(self, plotly_name="selectionrevision", parent_name="layout", **kwargs):
        super(SelectionrevisionValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            **kwargs,
        )
