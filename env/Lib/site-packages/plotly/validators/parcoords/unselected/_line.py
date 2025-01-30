import _plotly_utils.basevalidators


class LineValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(
        self, plotly_name="line", parent_name="parcoords.unselected", **kwargs
    ):
        super(LineValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Line"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            color
                Sets the base color of unselected lines. in
                connection with `unselected.line.opacity`.
            opacity
                Sets the opacity of unselected lines. The
                default "auto" decreases the opacity smoothly
                as the number of lines increases. Use 1 to
                achieve exact `unselected.line.color`.
""",
            ),
            **kwargs,
        )
