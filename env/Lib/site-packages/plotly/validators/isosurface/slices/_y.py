import _plotly_utils.basevalidators


class YValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="y", parent_name="isosurface.slices", **kwargs):
        super(YValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Y"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            fill
                Sets the fill ratio of the `slices`. The
                default fill value of the `slices` is 1 meaning
                that they are entirely shaded. On the other
                hand Applying a `fill` ratio less than one
                would allow the creation of openings parallel
                to the edges.
            locations
                Specifies the location(s) of slices on the
                axis. When not specified slices would be
                created for all points of the axis y except
                start and end.
            locationssrc
                Sets the source reference on Chart Studio Cloud
                for `locations`.
            show
                Determines whether or not slice planes about
                the y dimension are drawn.
""",
            ),
            **kwargs,
        )
