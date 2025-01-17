import _plotly_utils.basevalidators


class ImagValidator(_plotly_utils.basevalidators.DataArrayValidator):
    def __init__(self, plotly_name="imag", parent_name="scattersmith", **kwargs):
        super(ImagValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc+clearAxisTypes"),
            **kwargs,
        )
