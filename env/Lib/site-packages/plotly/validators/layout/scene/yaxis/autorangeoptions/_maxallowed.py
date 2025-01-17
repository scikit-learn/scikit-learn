import _plotly_utils.basevalidators


class MaxallowedValidator(_plotly_utils.basevalidators.AnyValidator):
    def __init__(
        self,
        plotly_name="maxallowed",
        parent_name="layout.scene.yaxis.autorangeoptions",
        **kwargs,
    ):
        super(MaxallowedValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            implied_edits=kwargs.pop("implied_edits", {}),
            **kwargs,
        )
