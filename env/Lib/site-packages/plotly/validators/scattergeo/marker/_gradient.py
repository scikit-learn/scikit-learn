import _plotly_utils.basevalidators


class GradientValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(
        self, plotly_name="gradient", parent_name="scattergeo.marker", **kwargs
    ):
        super(GradientValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Gradient"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            color
                Sets the final color of the gradient fill: the
                center color for radial, the right for
                horizontal, or the bottom for vertical.
            colorsrc
                Sets the source reference on Chart Studio Cloud
                for `color`.
            type
                Sets the type of gradient used to fill the
                markers
            typesrc
                Sets the source reference on Chart Studio Cloud
                for `type`.
""",
            ),
            **kwargs,
        )
