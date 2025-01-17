import _plotly_utils.basevalidators


class SceneValidator(_plotly_utils.basevalidators.SubplotidValidator):
    def __init__(self, plotly_name="scene", parent_name="mesh3d", **kwargs):
        super(SceneValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            dflt=kwargs.pop("dflt", "scene"),
            edit_type=kwargs.pop("edit_type", "calc+clearAxisTypes"),
            **kwargs,
        )
