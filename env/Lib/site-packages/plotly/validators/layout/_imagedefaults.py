import _plotly_utils.basevalidators


class ImagedefaultsValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="imagedefaults", parent_name="layout", **kwargs):
        super(ImagedefaultsValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Image"),
            data_docs=kwargs.pop(
                "data_docs",
                """
""",
            ),
            **kwargs,
        )
