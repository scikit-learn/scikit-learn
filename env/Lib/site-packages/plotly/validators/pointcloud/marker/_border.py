import _plotly_utils.basevalidators


class BorderValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="border", parent_name="pointcloud.marker", **kwargs):
        super(BorderValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Border"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            arearatio
                Specifies what fraction of the marker area is
                covered with the border.
            color
                Sets the stroke color. It accepts a specific
                color. If the color is not fully opaque and
                there are hundreds of thousands of points, it
                may cause slower zooming and panning.
""",
            ),
            **kwargs,
        )
