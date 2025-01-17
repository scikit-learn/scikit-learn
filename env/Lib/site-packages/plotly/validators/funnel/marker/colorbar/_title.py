import _plotly_utils.basevalidators


class TitleValidator(_plotly_utils.basevalidators.TitleValidator):
    def __init__(
        self, plotly_name="title", parent_name="funnel.marker.colorbar", **kwargs
    ):
        super(TitleValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Title"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            font
                Sets this color bar's title font. Note that the
                title's font used to be set by the now
                deprecated `titlefont` attribute.
            side
                Determines the location of color bar's title
                with respect to the color bar. Defaults to
                "top" when `orientation` if "v" and  defaults
                to "right" when `orientation` if "h". Note that
                the title's location used to be set by the now
                deprecated `titleside` attribute.
            text
                Sets the title of the color bar. Note that
                before the existence of `title.text`, the
                title's contents used to be defined as the
                `title` attribute itself. This behavior has
                been deprecated.
""",
            ),
            **kwargs,
        )
