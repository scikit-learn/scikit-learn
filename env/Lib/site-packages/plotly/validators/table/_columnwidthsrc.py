import _plotly_utils.basevalidators


class ColumnwidthsrcValidator(_plotly_utils.basevalidators.SrcValidator):
    def __init__(self, plotly_name="columnwidthsrc", parent_name="table", **kwargs):
        super(ColumnwidthsrcValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "none"),
            **kwargs,
        )
