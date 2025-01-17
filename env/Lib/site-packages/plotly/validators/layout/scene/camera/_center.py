import _plotly_utils.basevalidators


class CenterValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(
        self, plotly_name="center", parent_name="layout.scene.camera", **kwargs
    ):
        super(CenterValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Center"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            x

            y

            z

""",
            ),
            **kwargs,
        )
