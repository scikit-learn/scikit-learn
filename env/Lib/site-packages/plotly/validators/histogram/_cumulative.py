import _plotly_utils.basevalidators


class CumulativeValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="cumulative", parent_name="histogram", **kwargs):
        super(CumulativeValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Cumulative"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            currentbin
                Only applies if cumulative is enabled. Sets
                whether the current bin is included, excluded,
                or has half of its value included in the
                current cumulative value. "include" is the
                default for compatibility with various other
                tools, however it introduces a half-bin bias to
                the results. "exclude" makes the opposite half-
                bin bias, and "half" removes it.
            direction
                Only applies if cumulative is enabled. If
                "increasing" (default) we sum all prior bins,
                so the result increases from left to right. If
                "decreasing" we sum later bins so the result
                decreases from left to right.
            enabled
                If true, display the cumulative distribution by
                summing the binned values. Use the `direction`
                and `centralbin` attributes to tune the
                accumulation method. Note: in this mode, the
                "density" `histnorm` settings behave the same
                as their equivalents without "density": "" and
                "density" both rise to the number of data
                points, and "probability" and *probability
                density* both rise to the number of sample
                points.
""",
            ),
            **kwargs,
        )
