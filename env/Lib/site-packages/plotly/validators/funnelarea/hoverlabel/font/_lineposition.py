import _plotly_utils.basevalidators


class LinepositionValidator(_plotly_utils.basevalidators.FlaglistValidator):
    def __init__(
        self,
        plotly_name="lineposition",
        parent_name="funnelarea.hoverlabel.font",
        **kwargs,
    ):
        super(LinepositionValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            array_ok=kwargs.pop("array_ok", True),
            edit_type=kwargs.pop("edit_type", "none"),
            extras=kwargs.pop("extras", ["none"]),
            flags=kwargs.pop("flags", ["under", "over", "through"]),
            **kwargs,
        )
