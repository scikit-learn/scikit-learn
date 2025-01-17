import _plotly_utils.basevalidators


class SortpathsValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="sortpaths", parent_name="parcats", **kwargs):
        super(SortpathsValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            values=kwargs.pop("values", ["forward", "backward"]),
            **kwargs,
        )
