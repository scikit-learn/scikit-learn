import _plotly_utils.basevalidators


class SpaceframeValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="spaceframe", parent_name="volume", **kwargs):
        super(SpaceframeValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Spaceframe"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            fill
                Sets the fill ratio of the `spaceframe`
                elements. The default fill value is 1 meaning
                that they are entirely shaded. Applying a
                `fill` ratio less than one would allow the
                creation of openings parallel to the edges.
            show
                Displays/hides tetrahedron shapes between
                minimum and maximum iso-values. Often useful
                when either caps or surfaces are disabled or
                filled with values less than 1.
""",
            ),
            **kwargs,
        )
