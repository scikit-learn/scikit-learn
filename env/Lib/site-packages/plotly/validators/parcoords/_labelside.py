import _plotly_utils.basevalidators


class LabelsideValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="labelside", parent_name="parcoords", **kwargs):
        super(LabelsideValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            values=kwargs.pop("values", ["top", "bottom"]),
            **kwargs,
        )
