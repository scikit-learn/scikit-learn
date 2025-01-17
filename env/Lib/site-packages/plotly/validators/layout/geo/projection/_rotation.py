import _plotly_utils.basevalidators


class RotationValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(
        self, plotly_name="rotation", parent_name="layout.geo.projection", **kwargs
    ):
        super(RotationValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Rotation"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            lat
                Rotates the map along meridians (in degrees
                North).
            lon
                Rotates the map along parallels (in degrees
                East). Defaults to the center of the
                `lonaxis.range` values.
            roll
                Roll the map (in degrees) For example, a roll
                of 180 makes the map appear upside down.
""",
            ),
            **kwargs,
        )
