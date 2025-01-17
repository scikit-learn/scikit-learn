import _plotly_utils.basevalidators


class TitleValidator(_plotly_utils.basevalidators.TitleValidator):
    def __init__(self, plotly_name="title", parent_name="layout.xaxis", **kwargs):
        super(TitleValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Title"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            font
                Sets this axis' title font. Note that the
                title's font used to be customized by the now
                deprecated `titlefont` attribute.
            standoff
                Sets the standoff distance (in px) between the
                axis labels and the title text The default
                value is a function of the axis tick labels,
                the title `font.size` and the axis `linewidth`.
                Note that the axis title position is always
                constrained within the margins, so the actual
                standoff distance is always less than the set
                or default value. By setting `standoff` and
                turning on `automargin`, plotly.js will push
                the margins to fit the axis title at given
                standoff distance.
            text
                Sets the title of this axis. Note that before
                the existence of `title.text`, the title's
                contents used to be defined as the `title`
                attribute itself. This behavior has been
                deprecated.
""",
            ),
            **kwargs,
        )
