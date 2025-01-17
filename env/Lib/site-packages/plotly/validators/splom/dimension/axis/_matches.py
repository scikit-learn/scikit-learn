import _plotly_utils.basevalidators


class MatchesValidator(_plotly_utils.basevalidators.BooleanValidator):
    def __init__(
        self, plotly_name="matches", parent_name="splom.dimension.axis", **kwargs
    ):
        super(MatchesValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
