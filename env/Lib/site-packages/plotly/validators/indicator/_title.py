import _plotly_utils.basevalidators


class TitleValidator(_plotly_utils.basevalidators.TitleValidator):
    def __init__(self, plotly_name="title", parent_name="indicator", **kwargs):
        super(TitleValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Title"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            align
                Sets the horizontal alignment of the title. It
                defaults to `center` except for bullet charts
                for which it defaults to right.
            font
                Set the font used to display the title
            text
                Sets the title of this indicator.
""",
            ),
            **kwargs,
        )
