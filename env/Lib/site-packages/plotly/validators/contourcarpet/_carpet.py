import _plotly_utils.basevalidators


class CarpetValidator(_plotly_utils.basevalidators.StringValidator):
    def __init__(self, plotly_name="carpet", parent_name="contourcarpet", **kwargs):
        super(CarpetValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
