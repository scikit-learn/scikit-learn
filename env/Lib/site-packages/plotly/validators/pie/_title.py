import _plotly_utils.basevalidators


class TitleValidator(_plotly_utils.basevalidators.TitleValidator):
    def __init__(self, plotly_name="title", parent_name="pie", **kwargs):
        super(TitleValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Title"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            font
                Sets the font used for `title`. Note that the
                title's font used to be set by the now
                deprecated `titlefont` attribute.
            position
                Specifies the location of the `title`. Note
                that the title's position used to be set by the
                now deprecated `titleposition` attribute.
            text
                Sets the title of the chart. If it is empty, no
                title is displayed. Note that before the
                existence of `title.text`, the title's contents
                used to be defined as the `title` attribute
                itself. This behavior has been deprecated.
""",
            ),
            **kwargs,
        )
