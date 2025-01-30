import _plotly_utils.basevalidators


class Copy_ZstyleValidator(_plotly_utils.basevalidators.BooleanValidator):
    def __init__(
        self, plotly_name="copy_zstyle", parent_name="scatter3d.error_x", **kwargs
    ):
        super(Copy_ZstyleValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "calc"),
            **kwargs,
        )
