import _plotly_utils.basevalidators


class HoverlabelValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(
        self, plotly_name="hoverlabel", parent_name="layout.scene.annotation", **kwargs
    ):
        super(HoverlabelValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Hoverlabel"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            bgcolor
                Sets the background color of the hover label.
                By default uses the annotation's `bgcolor` made
                opaque, or white if it was transparent.
            bordercolor
                Sets the border color of the hover label. By
                default uses either dark grey or white, for
                maximum contrast with `hoverlabel.bgcolor`.
            font
                Sets the hover label text font. By default uses
                the global hover font and size, with color from
                `hoverlabel.bordercolor`.
""",
            ),
            **kwargs,
        )
