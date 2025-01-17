import _plotly_utils.basevalidators


class ContourValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="contour", parent_name="mesh3d", **kwargs):
        super(ContourValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Contour"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            color
                Sets the color of the contour lines.
            show
                Sets whether or not dynamic contours are shown
                on hover
            width
                Sets the width of the contour lines.
""",
            ),
            **kwargs,
        )
