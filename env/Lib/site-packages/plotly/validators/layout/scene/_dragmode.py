import _plotly_utils.basevalidators


class DragmodeValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="dragmode", parent_name="layout.scene", **kwargs):
        super(DragmodeValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            values=kwargs.pop("values", ["orbit", "turntable", "zoom", "pan", False]),
            **kwargs,
        )
