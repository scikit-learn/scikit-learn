import _plotly_utils.basevalidators


class PatternValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="pattern", parent_name="bar.marker", **kwargs):
        super(PatternValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Pattern"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            bgcolor
                When there is no colorscale sets the color of
                background pattern fill. Defaults to a
                `marker.color` background when `fillmode` is
                "overlay". Otherwise, defaults to a transparent
                background.
            bgcolorsrc
                Sets the source reference on Chart Studio Cloud
                for `bgcolor`.
            fgcolor
                When there is no colorscale sets the color of
                foreground pattern fill. Defaults to a
                `marker.color` background when `fillmode` is
                "replace". Otherwise, defaults to dark grey or
                white to increase contrast with the `bgcolor`.
            fgcolorsrc
                Sets the source reference on Chart Studio Cloud
                for `fgcolor`.
            fgopacity
                Sets the opacity of the foreground pattern
                fill. Defaults to a 0.5 when `fillmode` is
                "overlay". Otherwise, defaults to 1.
            fillmode
                Determines whether `marker.color` should be
                used as a default to `bgcolor` or a `fgcolor`.
            shape
                Sets the shape of the pattern fill. By default,
                no pattern is used for filling the area.
            shapesrc
                Sets the source reference on Chart Studio Cloud
                for `shape`.
            size
                Sets the size of unit squares of the pattern
                fill in pixels, which corresponds to the
                interval of repetition of the pattern.
            sizesrc
                Sets the source reference on Chart Studio Cloud
                for `size`.
            solidity
                Sets the solidity of the pattern fill. Solidity
                is roughly the fraction of the area filled by
                the pattern. Solidity of 0 shows only the
                background color without pattern and solidty of
                1 shows only the foreground color without
                pattern.
            soliditysrc
                Sets the source reference on Chart Studio Cloud
                for `solidity`.
""",
            ),
            **kwargs,
        )
