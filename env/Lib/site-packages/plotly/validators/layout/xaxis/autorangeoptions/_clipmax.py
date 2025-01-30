import _plotly_utils.basevalidators


class ClipmaxValidator(_plotly_utils.basevalidators.AnyValidator):
    def __init__(
        self,
        plotly_name="clipmax",
        parent_name="layout.xaxis.autorangeoptions",
        **kwargs,
    ):
        super(ClipmaxValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            implied_edits=kwargs.pop("implied_edits", {}),
            **kwargs,
        )
