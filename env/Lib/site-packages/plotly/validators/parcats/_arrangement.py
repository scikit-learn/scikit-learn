import _plotly_utils.basevalidators


class ArrangementValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="arrangement", parent_name="parcats", **kwargs):
        super(ArrangementValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            values=kwargs.pop("values", ["perpendicular", "freeform", "fixed"]),
            **kwargs,
        )
