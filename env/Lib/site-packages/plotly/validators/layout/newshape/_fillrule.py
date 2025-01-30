import _plotly_utils.basevalidators


class FillruleValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="fillrule", parent_name="layout.newshape", **kwargs):
        super(FillruleValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            values=kwargs.pop("values", ["evenodd", "nonzero"]),
            **kwargs,
        )
