import _plotly_utils.basevalidators


class AnglerefValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="angleref", parent_name="scatter.marker", **kwargs):
        super(AnglerefValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            anim=kwargs.pop("anim", False),
            edit_type=kwargs.pop("edit_type", "plot"),
            values=kwargs.pop("values", ["previous", "up"]),
            **kwargs,
        )
