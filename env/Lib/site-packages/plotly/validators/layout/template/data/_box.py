import _plotly_utils.basevalidators


class BoxValidator(_plotly_utils.basevalidators.CompoundArrayValidator):
    def __init__(self, plotly_name="box", parent_name="layout.template.data", **kwargs):
        super(BoxValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Box"),
            data_docs=kwargs.pop(
                "data_docs",
                """
""",
            ),
            **kwargs,
        )
