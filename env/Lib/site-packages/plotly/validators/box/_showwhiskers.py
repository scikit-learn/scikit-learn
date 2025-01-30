import _plotly_utils.basevalidators


class ShowwhiskersValidator(_plotly_utils.basevalidators.BooleanValidator):
    def __init__(self, plotly_name="showwhiskers", parent_name="box", **kwargs):
        super(ShowwhiskersValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
