import _plotly_utils.basevalidators


class MeanlineValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="meanline", parent_name="violin", **kwargs):
        super(MeanlineValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Meanline"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            color
                Sets the mean line color.
            visible
                Determines if a line corresponding to the
                sample's mean is shown inside the violins. If
                `box.visible` is turned on, the mean line is
                drawn inside the inner box. Otherwise, the mean
                line is drawn from one side of the violin to
                other.
            width
                Sets the mean line width.
""",
            ),
            **kwargs,
        )
