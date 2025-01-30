import _plotly_utils.basevalidators


class FontValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(
        self, plotly_name="font", parent_name="isosurface.hoverlabel", **kwargs
    ):
        super(FontValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Font"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            color

            colorsrc
                Sets the source reference on Chart Studio Cloud
                for `color`.
            family
                HTML font family - the typeface that will be
                applied by the web browser. The web browser
                will only be able to apply a font if it is
                available on the system which it operates.
                Provide multiple font families, separated by
                commas, to indicate the preference in which to
                apply fonts if they aren't available on the
                system. The Chart Studio Cloud (at
                https://chart-studio.plotly.com or on-premise)
                generates images on a server, where only a
                select number of fonts are installed and
                supported. These include "Arial", "Balto",
                "Courier New", "Droid Sans", "Droid Serif",
                "Droid Sans Mono", "Gravitas One", "Old
                Standard TT", "Open Sans", "Overpass", "PT Sans
                Narrow", "Raleway", "Times New Roman".
            familysrc
                Sets the source reference on Chart Studio Cloud
                for `family`.
            lineposition
                Sets the kind of decoration line(s) with text,
                such as an "under", "over" or "through" as well
                as combinations e.g. "under+over", etc.
            linepositionsrc
                Sets the source reference on Chart Studio Cloud
                for `lineposition`.
            shadow
                Sets the shape and color of the shadow behind
                text. "auto" places minimal shadow and applies
                contrast text font color. See
                https://developer.mozilla.org/en-
                US/docs/Web/CSS/text-shadow for additional
                options.
            shadowsrc
                Sets the source reference on Chart Studio Cloud
                for `shadow`.
            size

            sizesrc
                Sets the source reference on Chart Studio Cloud
                for `size`.
            style
                Sets whether a font should be styled with a
                normal or italic face from its family.
            stylesrc
                Sets the source reference on Chart Studio Cloud
                for `style`.
            textcase
                Sets capitalization of text. It can be used to
                make text appear in all-uppercase or all-
                lowercase, or with each word capitalized.
            textcasesrc
                Sets the source reference on Chart Studio Cloud
                for `textcase`.
            variant
                Sets the variant of the font.
            variantsrc
                Sets the source reference on Chart Studio Cloud
                for `variant`.
            weight
                Sets the weight (or boldness) of the font.
            weightsrc
                Sets the source reference on Chart Studio Cloud
                for `weight`.
""",
            ),
            **kwargs,
        )
