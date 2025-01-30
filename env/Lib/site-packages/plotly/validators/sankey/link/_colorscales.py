import _plotly_utils.basevalidators


class ColorscalesValidator(_plotly_utils.basevalidators.CompoundArrayValidator):
    def __init__(self, plotly_name="colorscales", parent_name="sankey.link", **kwargs):
        super(ColorscalesValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Colorscale"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            cmax
                Sets the upper bound of the color domain.
            cmin
                Sets the lower bound of the color domain.
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
            label
                The label of the links to color based on their
                concentration within a flow.
            name
                When used in a template, named items are
                created in the output figure in addition to any
                items the figure already has in this array. You
                can modify these items in the output figure by
                making your own item with `templateitemname`
                matching this `name` alongside your
                modifications (including `visible: false` or
                `enabled: false` to hide it). Has no effect
                outside of a template.
            templateitemname
                Used to refer to a named item in this array in
                the template. Named items from the template
                will be created even without a matching item in
                the input figure, but you can modify one by
                making an item with `templateitemname` matching
                its `name`, alongside your modifications
                (including `visible: false` or `enabled: false`
                to hide it). If there is no template or no
                matching item, this item will be hidden unless
                you explicitly show it with `visible: true`.
""",
            ),
            **kwargs,
        )
