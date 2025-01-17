import _plotly_utils.basevalidators


class MarginValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="margin", parent_name="layout", **kwargs):
        super(MarginValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Margin"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            autoexpand
                Turns on/off margin expansion computations.
                Legends, colorbars, updatemenus, sliders, axis
                rangeselector and rangeslider are allowed to
                push the margins by defaults.
            b
                Sets the bottom margin (in px).
            l
                Sets the left margin (in px).
            pad
                Sets the amount of padding (in px) between the
                plotting area and the axis lines
            r
                Sets the right margin (in px).
            t
                Sets the top margin (in px).
""",
            ),
            **kwargs,
        )
