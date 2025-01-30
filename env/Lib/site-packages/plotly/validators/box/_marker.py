import _plotly_utils.basevalidators


class MarkerValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="marker", parent_name="box", **kwargs):
        super(MarkerValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Marker"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            angle
                Sets the marker angle in respect to `angleref`.
            color
                Sets the marker color. It accepts either a
                specific color or an array of numbers that are
                mapped to the colorscale relative to the max
                and min values of the array or relative to
                `marker.cmin` and `marker.cmax` if set.
            line
                :class:`plotly.graph_objects.box.marker.Line`
                instance or dict with compatible properties
            opacity
                Sets the marker opacity.
            outliercolor
                Sets the color of the outlier sample points.
            size
                Sets the marker size (in px).
            symbol
                Sets the marker symbol type. Adding 100 is
                equivalent to appending "-open" to a symbol
                name. Adding 200 is equivalent to appending
                "-dot" to a symbol name. Adding 300 is
                equivalent to appending "-open-dot" or "dot-
                open" to a symbol name.
""",
            ),
            **kwargs,
        )
