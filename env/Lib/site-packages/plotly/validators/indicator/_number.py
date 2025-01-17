import _plotly_utils.basevalidators


class NumberValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="number", parent_name="indicator", **kwargs):
        super(NumberValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Number"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            font
                Set the font used to display main number
            prefix
                Sets a prefix appearing before the number.
            suffix
                Sets a suffix appearing next to the number.
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
