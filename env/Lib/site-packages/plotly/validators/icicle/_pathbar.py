import _plotly_utils.basevalidators


class PathbarValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="pathbar", parent_name="icicle", **kwargs):
        super(PathbarValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Pathbar"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            edgeshape
                Determines which shape is used for edges
                between `barpath` labels.
            side
                Determines on which side of the the treemap the
                `pathbar` should be presented.
            textfont
                Sets the font used inside `pathbar`.
            thickness
                Sets the thickness of `pathbar` (in px). If not
                specified the `pathbar.textfont.size` is used
                with 3 pixles extra padding on each side.
            visible
                Determines if the path bar is drawn i.e.
                outside the trace `domain` and with one pixel
                gap.
""",
            ),
            **kwargs,
        )
