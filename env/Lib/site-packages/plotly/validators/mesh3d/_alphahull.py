import _plotly_utils.basevalidators


class AlphahullValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(self, plotly_name="alphahull", parent_name="mesh3d", **kwargs):
        super(AlphahullValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
