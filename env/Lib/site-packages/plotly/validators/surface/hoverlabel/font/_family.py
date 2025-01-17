import _plotly_utils.basevalidators


class FamilyValidator(_plotly_utils.basevalidators.StringValidator):
    def __init__(
        self, plotly_name="family", parent_name="surface.hoverlabel.font", **kwargs
    ):
        super(FamilyValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            array_ok=kwargs.pop("array_ok", True),
            edit_type=kwargs.pop("edit_type", "none"),
            no_blank=kwargs.pop("no_blank", True),
            strict=kwargs.pop("strict", True),
            **kwargs,
        )
