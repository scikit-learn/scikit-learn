import _plotly_utils.basevalidators


class ZValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="z", parent_name="volume.caps", **kwargs):
        super(ZValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Z"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            fill
                Sets the fill ratio of the `caps`. The default
                fill value of the `caps` is 1 meaning that they
                are entirely shaded. On the other hand Applying
                a `fill` ratio less than one would allow the
                creation of openings parallel to the edges.
            show
                Sets the fill ratio of the `slices`. The
                default fill value of the z `slices` is 1
                meaning that they are entirely shaded. On the
                other hand Applying a `fill` ratio less than
                one would allow the creation of openings
                parallel to the edges.
""",
            ),
            **kwargs,
        )
