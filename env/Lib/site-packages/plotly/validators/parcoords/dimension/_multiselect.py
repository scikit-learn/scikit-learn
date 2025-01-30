import _plotly_utils.basevalidators


class MultiselectValidator(_plotly_utils.basevalidators.BooleanValidator):
    def __init__(
        self, plotly_name="multiselect", parent_name="parcoords.dimension", **kwargs
    ):
        super(MultiselectValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            **kwargs,
        )
