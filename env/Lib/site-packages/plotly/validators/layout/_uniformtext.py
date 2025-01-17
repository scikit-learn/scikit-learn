import _plotly_utils.basevalidators


class UniformtextValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="uniformtext", parent_name="layout", **kwargs):
        super(UniformtextValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Uniformtext"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            minsize
                Sets the minimum text size between traces of
                the same type.
            mode
                Determines how the font size for various text
                elements are uniformed between each trace type.
                If the computed text sizes were smaller than
                the minimum size defined by
                `uniformtext.minsize` using "hide" option hides
                the text; and using "show" option shows the
                text without further downscaling. Please note
                that if the size defined by `minsize` is
                greater than the font size defined by trace,
                then the `minsize` is used.
""",
            ),
            **kwargs,
        )
