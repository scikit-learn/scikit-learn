import _plotly_utils.basevalidators


class SurfaceValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="surface", parent_name="isosurface", **kwargs):
        super(SurfaceValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Surface"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            count
                Sets the number of iso-surfaces between minimum
                and maximum iso-values. By default this value
                is 2 meaning that only minimum and maximum
                surfaces would be drawn.
            fill
                Sets the fill ratio of the iso-surface. The
                default fill value of the surface is 1 meaning
                that they are entirely shaded. On the other
                hand Applying a `fill` ratio less than one
                would allow the creation of openings parallel
                to the edges.
            pattern
                Sets the surface pattern of the iso-surface 3-D
                sections. The default pattern of the surface is
                `all` meaning that the rest of surface elements
                would be shaded. The check options (either 1 or
                2) could be used to draw half of the squares on
                the surface. Using various combinations of
                capital `A`, `B`, `C`, `D` and `E` may also be
                used to reduce the number of triangles on the
                iso-surfaces and creating other patterns of
                interest.
            show
                Hides/displays surfaces between minimum and
                maximum iso-values.
""",
            ),
            **kwargs,
        )
