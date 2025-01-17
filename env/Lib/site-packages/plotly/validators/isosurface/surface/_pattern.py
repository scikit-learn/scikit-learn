import _plotly_utils.basevalidators


class PatternValidator(_plotly_utils.basevalidators.FlaglistValidator):
    def __init__(
        self, plotly_name="pattern", parent_name="isosurface.surface", **kwargs
    ):
        super(PatternValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            extras=kwargs.pop("extras", ["all", "odd", "even"]),
            flags=kwargs.pop("flags", ["A", "B", "C", "D", "E"]),
            **kwargs,
        )
