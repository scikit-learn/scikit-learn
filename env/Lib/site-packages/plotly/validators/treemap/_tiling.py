import _plotly_utils.basevalidators


class TilingValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="tiling", parent_name="treemap", **kwargs):
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
            packing
                Determines d3 treemap solver. For more info
                please refer to
                https://github.com/d3/d3-hierarchy#treemap-
                tiling
            pad
                Sets the inner padding (in px).
            squarifyratio
                When using "squarify" `packing` algorithm,
                according to https://github.com/d3/d3-
                hierarchy/blob/v3.1.1/README.md#squarify_ratio
                this option specifies the desired aspect ratio
                of the generated rectangles. The ratio must be
                specified as a number greater than or equal to
                one. Note that the orientation of the generated
                rectangles (tall or wide) is not implied by the
                ratio; for example, a ratio of two will attempt
                to produce a mixture of rectangles whose
                width:height ratio is either 2:1 or 1:2. When
                using "squarify", unlike d3 which uses the
                Golden Ratio i.e. 1.618034, Plotly applies 1 to
                increase squares in treemap layouts.
""",
            ),
            **kwargs,
        )
