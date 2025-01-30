import _plotly_utils.basevalidators


class CircleValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(
        self, plotly_name="circle", parent_name="layout.mapbox.layer", **kwargs
    ):
        super(CircleValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Circle"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            radius
                Sets the circle radius
                (mapbox.layer.paint.circle-radius). Has an
                effect only when `type` is set to "circle".
""",
            ),
            **kwargs,
        )
