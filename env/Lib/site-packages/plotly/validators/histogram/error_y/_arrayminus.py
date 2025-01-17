import _plotly_utils.basevalidators


class ArrayminusValidator(_plotly_utils.basevalidators.DataArrayValidator):
    def __init__(
        self, plotly_name="arrayminus", parent_name="histogram.error_y", **kwargs
    ):
        super(ArrayminusValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
