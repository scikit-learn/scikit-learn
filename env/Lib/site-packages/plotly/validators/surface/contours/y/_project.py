import _plotly_utils.basevalidators


class ProjectValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(
        self, plotly_name="project", parent_name="surface.contours.y", **kwargs
    ):
        super(ProjectValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Project"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            x
                Determines whether or not these contour lines
                are projected on the x plane. If `highlight` is
                set to True (the default), the projected lines
                are shown on hover. If `show` is set to True,
                the projected lines are shown in permanence.
            y
                Determines whether or not these contour lines
                are projected on the y plane. If `highlight` is
                set to True (the default), the projected lines
                are shown on hover. If `show` is set to True,
                the projected lines are shown in permanence.
            z
                Determines whether or not these contour lines
                are projected on the z plane. If `highlight` is
                set to True (the default), the projected lines
                are shown on hover. If `show` is set to True,
                the projected lines are shown in permanence.
""",
            ),
            **kwargs,
        )
