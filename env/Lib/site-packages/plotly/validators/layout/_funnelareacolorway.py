import _plotly_utils.basevalidators


class FunnelareacolorwayValidator(_plotly_utils.basevalidators.ColorlistValidator):
    def __init__(
        self, plotly_name="funnelareacolorway", parent_name="layout", **kwargs
    ):
        super(FunnelareacolorwayValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
