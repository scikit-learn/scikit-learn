import _plotly_utils.basevalidators


class GroupclickValidator(_plotly_utils.basevalidators.EnumeratedValidator):
    def __init__(self, plotly_name="groupclick", parent_name="layout.legend", **kwargs):
        super(GroupclickValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "legend"),
            values=kwargs.pop("values", ["toggleitem", "togglegroup"]),
            **kwargs,
        )
