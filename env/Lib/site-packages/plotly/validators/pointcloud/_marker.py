import _plotly_utils.basevalidators


class MarkerValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="marker", parent_name="pointcloud", **kwargs):
        super(MarkerValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Marker"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            blend
                Determines if colors are blended together for a
                translucency effect in case `opacity` is
                specified as a value less then `1`. Setting
                `blend` to `true` reduces zoom/pan speed if
                used with large numbers of points.
            border
                :class:`plotly.graph_objects.pointcloud.marker.
                Border` instance or dict with compatible
                properties
            color
                Sets the marker fill color. It accepts a
                specific color. If the color is not fully
                opaque and there are hundreds of thousands of
                points, it may cause slower zooming and
                panning.
            opacity
                Sets the marker opacity. The default value is
                `1` (fully opaque). If the markers are not
                fully opaque and there are hundreds of
                thousands of points, it may cause slower
                zooming and panning. Opacity fades the color
                even if `blend` is left on `false` even if
                there is no translucency effect in that case.
            sizemax
                Sets the maximum size (in px) of the rendered
                marker points. Effective when the `pointcloud`
                shows only few points.
            sizemin
                Sets the minimum size (in px) of the rendered
                marker points, effective when the `pointcloud`
                shows a million or more points.
""",
            ),
            **kwargs,
        )
