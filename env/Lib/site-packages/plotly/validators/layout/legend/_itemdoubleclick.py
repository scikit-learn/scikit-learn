import _plotly_utils.basevalidators


class ItemdoubleclickValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(
        self, plotly_name="itemdoubleclick", parent_name="layout.legend", **kwargs
    ):
        super(ItemdoubleclickValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "legend"),
            values=kwargs.pop("values", ["toggle", "toggleothers", False]),
            **kwargs,
        )
