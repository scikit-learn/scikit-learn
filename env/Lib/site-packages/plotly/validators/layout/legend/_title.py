import _plotly_utils.basevalidators


class TitleValidator(_plotly_utils.basevalidators.TitleValidator):
    def __init__(self, plotly_name="title", parent_name="layout.legend", **kwargs):
        super(TitleValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Title"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            font
                Sets this legend's title font. Defaults to
                `legend.font` with its size increased about
                20%.
            side
                Determines the location of legend's title with
                respect to the legend items. Defaulted to "top"
                with `orientation` is "h". Defaulted to "left"
                with `orientation` is "v". The *top left*
                options could be used to expand top center and
                top right are for horizontal alignment legend
                area in both x and y sides.
            text
                Sets the title of the legend.
""",
            ),
            **kwargs,
        )
