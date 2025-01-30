import _plotly_utils.basevalidators


class ColoraxisValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="coloraxis", parent_name="layout", **kwargs):
        super(ColoraxisValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Coloraxis"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            autocolorscale
                Determines whether the colorscale is a default
                palette (`autocolorscale: true`) or the palette
                determined by `colorscale`. In case
                `colorscale` is unspecified or `autocolorscale`
                is true, the default palette will be chosen
                according to whether numbers in the `color`
                array are all positive, all negative or mixed.
            cauto
                Determines whether or not the color domain is
                computed with respect to the input data (here
                corresponding trace color array(s)) or the
                bounds set in `cmin` and `cmax` Defaults to
                `false` when `cmin` and `cmax` are set by the
                user.
            cmax
                Sets the upper bound of the color domain. Value
                should have the same units as corresponding
                trace color array(s) and if set, `cmin` must be
                set as well.
            cmid
                Sets the mid-point of the color domain by
                scaling `cmin` and/or `cmax` to be equidistant
                to this point. Value should have the same units
                as corresponding trace color array(s). Has no
                effect when `cauto` is `false`.
            cmin
                Sets the lower bound of the color domain. Value
                should have the same units as corresponding
                trace color array(s) and if set, `cmax` must be
                set as well.
            colorbar
                :class:`plotly.graph_objects.layout.coloraxis.C
                olorBar` instance or dict with compatible
                properties
            colorscale
                Sets the colorscale. The colorscale must be an
                array containing arrays mapping a normalized
                value to an rgb, rgba, hex, hsl, hsv, or named
                color string. At minimum, a mapping for the
                lowest (0) and highest (1) values are required.
                For example, `[[0, 'rgb(0,0,255)'], [1,
                'rgb(255,0,0)']]`. To control the bounds of the
                colorscale in color space, use `cmin` and
                `cmax`. Alternatively, `colorscale` may be a
                palette name string of the following list: Blac
                kbody,Bluered,Blues,Cividis,Earth,Electric,Gree
                ns,Greys,Hot,Jet,Picnic,Portland,Rainbow,RdBu,R
                eds,Viridis,YlGnBu,YlOrRd.
            reversescale
                Reverses the color mapping if true. If true,
                `cmin` will correspond to the last color in the
                array and `cmax` will correspond to the first
                color.
            showscale
                Determines whether or not a colorbar is
                displayed for this trace.
""",
            ),
            **kwargs,
        )
