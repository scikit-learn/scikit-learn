import _plotly_utils.basevalidators


class Arraytick0Validator(_plotly_utils.basevalidators.IntegerValidator):
    def __init__(self, plotly_name="arraytick0", parent_name="carpet.baxis", **kwargs):
        super(Arraytick0Validator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            min=kwargs.pop("min", 0),
            **kwargs,
        )
