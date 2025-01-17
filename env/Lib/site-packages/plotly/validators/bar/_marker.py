import _plotly_utils.basevalidators


class MarkerValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="marker", parent_name="bar", **kwargs):
        super(MarkerValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Marker"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            autocolorscale
                Determines whether the colorscale is a default
                palette (`autocolorscale: true`) or the palette
                determined by `marker.colorscale`. Has an
                effect only if in `marker.color` is set to a
                numerical array. In case `colorscale` is
                unspecified or `autocolorscale` is true, the
                default palette will be chosen according to
                whether numbers in the `color` array are all
                positive, all negative or mixed.
            cauto
                Determines whether or not the color domain is
                computed with respect to the input data (here
                in `marker.color`) or the bounds set in
                `marker.cmin` and `marker.cmax` Has an effect
                only if in `marker.color` is set to a numerical
                array. Defaults to `false` when `marker.cmin`
                and `marker.cmax` are set by the user.
            cmax
                Sets the upper bound of the color domain. Has
                an effect only if in `marker.color` is set to a
                numerical array. Value should have the same
                units as in `marker.color` and if set,
                `marker.cmin` must be set as well.
            cmid
                Sets the mid-point of the color domain by
                scaling `marker.cmin` and/or `marker.cmax` to
                be equidistant to this point. Has an effect
                only if in `marker.color` is set to a numerical
                array. Value should have the same units as in
                `marker.color`. Has no effect when
                `marker.cauto` is `false`.
            cmin
                Sets the lower bound of the color domain. Has
                an effect only if in `marker.color` is set to a
                numerical array. Value should have the same
                units as in `marker.color` and if set,
                `marker.cmax` must be set as well.
            color
                Sets the marker color. It accepts either a
                specific color or an array of numbers that are
                mapped to the colorscale relative to the max
                and min values of the array or relative to
                `marker.cmin` and `marker.cmax` if set.
            coloraxis
                Sets a reference to a shared color axis.
                References to these shared color axes are
                "coloraxis", "coloraxis2", "coloraxis3", etc.
                Settings for these shared color axes are set in
                the layout, under `layout.coloraxis`,
                `layout.coloraxis2`, etc. Note that multiple
                color scales can be linked to the same color
                axis.
            colorbar
                :class:`plotly.graph_objects.bar.marker.ColorBa
                r` instance or dict with compatible properties
            colorscale
                Sets the colorscale. Has an effect only if in
                `marker.color` is set to a numerical array. The
                colorscale must be an array containing arrays
                mapping a normalized value to an rgb, rgba,
                hex, hsl, hsv, or named color string. At
                minimum, a mapping for the lowest (0) and
                highest (1) values are required. For example,
                `[[0, 'rgb(0,0,255)'], [1, 'rgb(255,0,0)']]`.
                To control the bounds of the colorscale in
                color space, use `marker.cmin` and
                `marker.cmax`. Alternatively, `colorscale` may
                be a palette name string of the following list:
                Blackbody,Bluered,Blues,Cividis,Earth,Electric,
                Greens,Greys,Hot,Jet,Picnic,Portland,Rainbow,Rd
                Bu,Reds,Viridis,YlGnBu,YlOrRd.
            colorsrc
                Sets the source reference on Chart Studio Cloud
                for `color`.
            cornerradius
                Sets the rounding of corners. May be an integer
                number of pixels, or a percentage of bar width
                (as a string ending in %). Defaults to
                `layout.barcornerradius`. In stack or relative
                barmode, the first trace to set cornerradius is
                used for the whole stack.
            line
                :class:`plotly.graph_objects.bar.marker.Line`
                instance or dict with compatible properties
            opacity
                Sets the opacity of the bars.
            opacitysrc
                Sets the source reference on Chart Studio Cloud
                for `opacity`.
            pattern
                Sets the pattern within the marker.
            reversescale
                Reverses the color mapping if true. Has an
                effect only if in `marker.color` is set to a
                numerical array. If true, `marker.cmin` will
                correspond to the last color in the array and
                `marker.cmax` will correspond to the first
                color.
            showscale
                Determines whether or not a colorbar is
                displayed for this trace. Has an effect only if
                in `marker.color` is set to a numerical array.
""",
            ),
            **kwargs,
        )
