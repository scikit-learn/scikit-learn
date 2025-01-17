import _plotly_utils.basevalidators


class BoxValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="box", parent_name="violin", **kwargs):
        super(BoxValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Box"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            fillcolor
                Sets the inner box plot fill color.
            line
                :class:`plotly.graph_objects.violin.box.Line`
                instance or dict with compatible properties
            visible
                Determines if an miniature box plot is drawn
                inside the violins.
            width
                Sets the width of the inner box plots relative
                to the violins' width. For example, with 1, the
                inner box plots are as wide as the violins.
""",
            ),
            **kwargs,
        )
