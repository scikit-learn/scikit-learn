import _plotly_utils.basevalidators


class SubplotValidator(_plotly_utils.basevalidators.SubplotidValidator):
    def __init__(self, plotly_name="subplot", parent_name="choroplethmap", **kwargs):
        super(SubplotValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            dflt=kwargs.pop("dflt", "map"),
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
