import _plotly_utils.basevalidators


class FillgradientValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="fillgradient", parent_name="scatter", **kwargs):
        super(FillgradientValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Fillgradient"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            colorscale
                Sets the fill gradient colors as a color scale.
                The color scale is interpreted as a gradient
                applied in the direction specified by
                "orientation", from the lowest to the highest
                value of the scatter plot along that axis, or
                from the center to the most distant point from
                it, if orientation is "radial".
            start
                Sets the gradient start value. It is given as
                the absolute position on the axis determined by
                the orientiation. E.g., if orientation is
                "horizontal", the gradient will be horizontal
                and start from the x-position given by start.
                If omitted, the gradient starts at the lowest
                value of the trace along the respective axis.
                Ignored if orientation is "radial".
            stop
                Sets the gradient end value. It is given as the
                absolute position on the axis determined by the
                orientiation. E.g., if orientation is
                "horizontal", the gradient will be horizontal
                and end at the x-position given by end. If
                omitted, the gradient ends at the highest value
                of the trace along the respective axis. Ignored
                if orientation is "radial".
            type
                Sets the type/orientation of the color gradient
                for the fill. Defaults to "none".
""",
            ),
            **kwargs,
        )
