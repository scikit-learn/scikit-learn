import _plotly_utils.basevalidators


class AutotypenumbersValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(
        self, plotly_name="autotypenumbers", parent_name="layout.scene.yaxis", **kwargs
    ):
        super(AutotypenumbersValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            values=kwargs.pop("values", ["convert types", "strict"]),
            **kwargs,
        )
