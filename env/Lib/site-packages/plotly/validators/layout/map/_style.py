import _plotly_utils.basevalidators


class StyleValidator(_plotly_utils.basevalidators.AnyValidator):
    def __init__(self, plotly_name="style", parent_name="layout.map", **kwargs):
        super(StyleValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            values=kwargs.pop(
                "values",
                [
                    "basic",
                    "carto-darkmatter",
                    "carto-darkmatter-nolabels",
                    "carto-positron",
                    "carto-positron-nolabels",
                    "carto-voyager",
                    "carto-voyager-nolabels",
                    "dark",
                    "light",
                    "open-street-map",
                    "outdoors",
                    "satellite",
                    "satellite-streets",
                    "streets",
                    "white-bg",
                ],
            ),
            **kwargs,
        )
