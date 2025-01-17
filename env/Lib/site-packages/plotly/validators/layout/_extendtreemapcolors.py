import _plotly_utils.basevalidators


class ExtendtreemapcolorsValidator(_plotly_utils.basevalidators.BooleanValidator):
    def __init__(
        self, plotly_name="extendtreemapcolors", parent_name="layout", **kwargs
    ):
        super(ExtendtreemapcolorsValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
