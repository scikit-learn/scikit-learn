import _plotly_utils.basevalidators


class AnglerefValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(
        self, plotly_name="angleref", parent_name="scattersmith.marker", **kwargs
    ):
        super(AnglerefValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            values=kwargs.pop("values", ["previous", "up"]),
            **kwargs,
        )
