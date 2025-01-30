import _plotly_utils.basevalidators


class LineValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="line", parent_name="layout.mapbox.layer", **kwargs):
        super(LineValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Line"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            dash
                Sets the length of dashes and gaps
                (mapbox.layer.paint.line-dasharray). Has an
                effect only when `type` is set to "line".
            dashsrc
                Sets the source reference on Chart Studio Cloud
                for `dash`.
            width
                Sets the line width (mapbox.layer.paint.line-
                width). Has an effect only when `type` is set
                to "line".
""",
            ),
            **kwargs,
        )
