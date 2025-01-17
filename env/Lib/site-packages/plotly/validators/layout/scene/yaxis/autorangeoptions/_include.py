import _plotly_utils.basevalidators


class IncludeValidator(_plotly_utils.basevalidators.AnyValidator):
    def __init__(
        self,
        plotly_name="include",
        parent_name="layout.scene.yaxis.autorangeoptions",
        **kwargs,
    ):
        super(IncludeValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            array_ok=kwargs.pop("array_ok", True),
            edit_type=kwargs.pop("edit_type", "plot"),
            implied_edits=kwargs.pop("implied_edits", {}),
            **kwargs,
        )
