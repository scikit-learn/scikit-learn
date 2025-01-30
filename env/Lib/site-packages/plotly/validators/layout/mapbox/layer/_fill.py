import _plotly_utils.basevalidators


class FillValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="fill", parent_name="layout.mapbox.layer", **kwargs):
        super(FillValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Fill"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            outlinecolor
                Sets the fill outline color
                (mapbox.layer.paint.fill-outline-color). Has an
                effect only when `type` is set to "fill".
""",
            ),
            **kwargs,
        )
