import _plotly_utils.basevalidators


class TextinfoValidator(_plotly_utils.basevalidators.FlaglistValidator):
    def __init__(self, plotly_name="textinfo", parent_name="sunburst", **kwargs):
        super(TextinfoValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            extras=kwargs.pop("extras", ["none"]),
            flags=kwargs.pop(
                "flags",
                [
                    "label",
                    "text",
                    "value",
                    "current path",
                    "percent root",
                    "percent entry",
                    "percent parent",
                ],
            ),
            **kwargs,
        )
