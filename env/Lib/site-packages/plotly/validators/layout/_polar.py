import _plotly_utils.basevalidators


class PolarValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="polar", parent_name="layout", **kwargs):
        super(PolarValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Polar"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            angularaxis
                :class:`plotly.graph_objects.layout.polar.Angul
                arAxis` instance or dict with compatible
                properties
            bargap
                Sets the gap between bars of adjacent location
                coordinates. Values are unitless, they
                represent fractions of the minimum difference
                in bar positions in the data.
            barmode
                Determines how bars at the same location
                coordinate are displayed on the graph. With
                "stack", the bars are stacked on top of one
                another With "overlay", the bars are plotted
                over one another, you might need to reduce
                "opacity" to see multiple bars.
            bgcolor
                Set the background color of the subplot
            domain
                :class:`plotly.graph_objects.layout.polar.Domai
                n` instance or dict with compatible properties
            gridshape
                Determines if the radial axis grid lines and
                angular axis line are drawn as "circular"
                sectors or as "linear" (polygon) sectors. Has
                an effect only when the angular axis has `type`
                "category". Note that `radialaxis.angle` is
                snapped to the angle of the closest vertex when
                `gridshape` is "circular" (so that radial axis
                scale is the same as the data scale).
            hole
                Sets the fraction of the radius to cut out of
                the polar subplot.
            radialaxis
                :class:`plotly.graph_objects.layout.polar.Radia
                lAxis` instance or dict with compatible
                properties
            sector
                Sets angular span of this polar subplot with
                two angles (in degrees). Sector are assumed to
                be spanned in the counterclockwise direction
                with 0 corresponding to rightmost limit of the
                polar subplot.
            uirevision
                Controls persistence of user-driven changes in
                axis attributes, if not overridden in the
                individual axes. Defaults to
                `layout.uirevision`.
""",
            ),
            **kwargs,
        )
