import _plotly_utils.basevalidators


class FillValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="fill", parent_name="table.header", **kwargs):
        super(FillValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Fill"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            color
                Sets the cell fill color. It accepts either a
                specific color or an array of colors or a 2D
                array of colors.
            colorsrc
                Sets the source reference on Chart Studio Cloud
                for `color`.
""",
            ),
            **kwargs,
        )
