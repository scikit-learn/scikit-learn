import _plotly_utils.basevalidators


class TilingValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="tiling", parent_name="icicle", **kwargs):
        super(TilingValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Tiling"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            flip
                Determines if the positions obtained from
                solver are flipped on each axis.
            orientation
                When set in conjunction with `tiling.flip`,
                determines on which side the root nodes are
                drawn in the chart. If `tiling.orientation` is
                "v" and `tiling.flip` is "", the root nodes
                appear at the top. If `tiling.orientation` is
                "v" and `tiling.flip` is "y", the root nodes
                appear at the bottom. If `tiling.orientation`
                is "h" and `tiling.flip` is "", the root nodes
                appear at the left. If `tiling.orientation` is
                "h" and `tiling.flip` is "x", the root nodes
                appear at the right.
            pad
                Sets the inner padding (in px).
""",
            ),
            **kwargs,
        )
