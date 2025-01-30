import _plotly_utils.basevalidators


class ValuessrcValidator(_plotly_utils.basevalidators.SrcValidator):
    def __init__(
        self, plotly_name="valuessrc", parent_name="parcats.dimension", **kwargs
    ):
        super(ValuessrcValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            **kwargs,
        )
