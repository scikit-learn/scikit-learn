import _plotly_utils.basevalidators


class PadValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="pad", parent_name="treemap.marker", **kwargs):
        super(PadValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Pad"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            b
                Sets the padding form the bottom (in px).
            l
                Sets the padding form the left (in px).
            r
                Sets the padding form the right (in px).
            t
                Sets the padding form the top (in px).
""",
            ),
            **kwargs,
        )
