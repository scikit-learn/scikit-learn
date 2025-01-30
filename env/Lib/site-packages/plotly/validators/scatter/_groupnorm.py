import _plotly_utils.basevalidators


class GroupnormValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="groupnorm", parent_name="scatter", **kwargs):
        super(GroupnormValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            values=kwargs.pop("values", ["", "fraction", "percent"]),
            **kwargs,
        )
