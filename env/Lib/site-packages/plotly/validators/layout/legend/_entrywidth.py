import _plotly_utils.basevalidators


class EntrywidthValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="entrywidth", parent_name="layout.legend", **kwargs):
        super(EntrywidthValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "legend"),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
