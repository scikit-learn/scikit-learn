import _plotly_utils.basevalidators


class ShapedefaultsValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="shapedefaults", parent_name="layout", **kwargs):
        super(ShapedefaultsValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Shape"),
            data_docs=kwargs.pop(
                "data_docs",
                """
""",
            ),
            **kwargs,
        )
