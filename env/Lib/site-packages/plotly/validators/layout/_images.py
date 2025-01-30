import _plotly_utils.basevalidators


class ImagesValidator(_plotly_utils.basevalidators.CompoundArrayValidator):
    def __init__(self, plotly_name="images", parent_name="layout", **kwargs):
        super(ImagesValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Image"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            layer
                Specifies whether images are drawn below or
                above traces. When `xref` and `yref` are both
                set to `paper`, image is drawn below the entire
                plot area.
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
            opacity
                Sets the opacity of the image.
            sizex
                Sets the image container size horizontally. The
                image will be sized based on the `position`
                value. When `xref` is set to `paper`, units are
                sized relative to the plot width. When `xref`
                ends with ` domain`, units are sized relative
                to the axis width.
            sizey
                Sets the image container size vertically. The
                image will be sized based on the `position`
                value. When `yref` is set to `paper`, units are
                sized relative to the plot height. When `yref`
                ends with ` domain`, units are sized relative
                to the axis height.
            sizing
                Specifies which dimension of the image to
                constrain.
            source
                Specifies the URL of the image to be used. The
                URL must be accessible from the domain where
                the plot code is run, and can be either
                relative or absolute.
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
            visible
                Determines whether or not this image is
                visible.
            x
                Sets the image's x position. When `xref` is set
                to `paper`, units are sized relative to the
                plot height. See `xref` for more info
            xanchor
                Sets the anchor for the x position
            xref
                Sets the images's x coordinate axis. If set to
                a x axis id (e.g. "x" or "x2"), the `x`
                position refers to a x coordinate. If set to
                "paper", the `x` position refers to the
                distance from the left of the plotting area in
                normalized coordinates where 0 (1) corresponds
                to the left (right). If set to a x axis ID
                followed by "domain" (separated by a space),
                the position behaves like for "paper", but
                refers to the distance in fractions of the
                domain length from the left of the domain of
                that axis: e.g., *x2 domain* refers to the
                domain of the second x  axis and a x position
                of 0.5 refers to the point between the left and
                the right of the domain of the second x axis.
            y
                Sets the image's y position. When `yref` is set
                to `paper`, units are sized relative to the
                plot height. See `yref` for more info
            yanchor
                Sets the anchor for the y position.
            yref
                Sets the images's y coordinate axis. If set to
                a y axis id (e.g. "y" or "y2"), the `y`
                position refers to a y coordinate. If set to
                "paper", the `y` position refers to the
                distance from the bottom of the plotting area
                in normalized coordinates where 0 (1)
                corresponds to the bottom (top). If set to a y
                axis ID followed by "domain" (separated by a
                space), the position behaves like for "paper",
                but refers to the distance in fractions of the
                domain length from the bottom of the domain of
                that axis: e.g., *y2 domain* refers to the
                domain of the second y  axis and a y position
                of 0.5 refers to the point between the bottom
                and the top of the domain of the second y axis.
""",
            ),
            **kwargs,
        )
