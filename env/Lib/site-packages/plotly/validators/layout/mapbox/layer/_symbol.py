import _plotly_utils.basevalidators


class SymbolValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(
        self, plotly_name="symbol", parent_name="layout.mapbox.layer", **kwargs
    ):
        super(SymbolValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Symbol"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            icon
                Sets the symbol icon image
                (mapbox.layer.layout.icon-image). Full list:
                https://www.mapbox.com/maki-icons/
            iconsize
                Sets the symbol icon size
                (mapbox.layer.layout.icon-size). Has an effect
                only when `type` is set to "symbol".
            placement
                Sets the symbol and/or text placement
                (mapbox.layer.layout.symbol-placement). If
                `placement` is "point", the label is placed
                where the geometry is located If `placement` is
                "line", the label is placed along the line of
                the geometry If `placement` is "line-center",
                the label is placed on the center of the
                geometry
            text
                Sets the symbol text (mapbox.layer.layout.text-
                field).
            textfont
                Sets the icon text font
                (color=mapbox.layer.paint.text-color,
                size=mapbox.layer.layout.text-size). Has an
                effect only when `type` is set to "symbol".
            textposition
                Sets the positions of the `text` elements with
                respects to the (x,y) coordinates.
""",
            ),
            **kwargs,
        )
