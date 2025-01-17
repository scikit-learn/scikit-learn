import _plotly_utils.basevalidators


class MirrorValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(
        self, plotly_name="mirror", parent_name="layout.scene.xaxis", **kwargs
    ):
        super(MirrorValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            values=kwargs.pop("values", [True, "ticks", False, "all", "allticks"]),
            **kwargs,
        )
