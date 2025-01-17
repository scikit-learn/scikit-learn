import _plotly_utils.basevalidators


class CurrentvalueValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(
        self, plotly_name="currentvalue", parent_name="layout.slider", **kwargs
    ):
        super(CurrentvalueValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Currentvalue"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            font
                Sets the font of the current value label text.
            offset
                The amount of space, in pixels, between the
                current value label and the slider.
            prefix
                When currentvalue.visible is true, this sets
                the prefix of the label.
            suffix
                When currentvalue.visible is true, this sets
                the suffix of the label.
            visible
                Shows the currently-selected value above the
                slider.
            xanchor
                The alignment of the value readout relative to
                the length of the slider.
""",
            ),
            **kwargs,
        )
