import _plotly_utils.basevalidators


class DeltaValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="delta", parent_name="indicator", **kwargs):
        super(DeltaValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Delta"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            decreasing
                :class:`plotly.graph_objects.indicator.delta.De
                creasing` instance or dict with compatible
                properties
            font
                Set the font used to display the delta
            increasing
                :class:`plotly.graph_objects.indicator.delta.In
                creasing` instance or dict with compatible
                properties
            position
                Sets the position of delta with respect to the
                number.
            prefix
                Sets a prefix appearing before the delta.
            reference
                Sets the reference value to compute the delta.
                By default, it is set to the current value.
            relative
                Show relative change
            suffix
                Sets a suffix appearing next to the delta.
            valueformat
                Sets the value formatting rule using d3
                formatting mini-languages which are very
                similar to those in Python. For numbers, see: h
                ttps://github.com/d3/d3-format/tree/v1.4.5#d3-
                format.
""",
            ),
            **kwargs,
        )
