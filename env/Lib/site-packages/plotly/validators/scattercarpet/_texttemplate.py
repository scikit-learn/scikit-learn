import _plotly_utils.basevalidators


class TexttemplateValidator(_plotly_utils.basevalidators.StringValidator):
    def __init__(
        self, plotly_name="texttemplate", parent_name="scattercarpet", **kwargs
    ):
        super(TexttemplateValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            array_ok=kwargs.pop("array_ok", True),
            edit_type=kwargs.pop("edit_type", "plot"),
            **kwargs,
        )
