import _plotly_utils.basevalidators


class PackingValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="packing", parent_name="treemap.tiling", **kwargs):
        super(PackingValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            values=kwargs.pop(
                "values",
                ["squarify", "binary", "dice", "slice", "slice-dice", "dice-slice"],
            ),
            **kwargs,
        )
