import _plotly_utils.basevalidators


class PadValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="pad", parent_name="layout.updatemenu", **kwargs):
        super(PadValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Pad"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            b
                The amount of padding (in px) along the bottom
                of the component.
            l
                The amount of padding (in px) on the left side
                of the component.
            r
                The amount of padding (in px) on the right side
                of the component.
            t
                The amount of padding (in px) along the top of
                the component.
""",
            ),
            **kwargs,
        )
