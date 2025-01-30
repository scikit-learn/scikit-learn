import _plotly_utils.basevalidators


class TemplateValidator(_plotly_utils.basevalidators.BaseTemplateValidator):
    def __init__(self, plotly_name="template", parent_name="layout", **kwargs):
        super(TemplateValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Template"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            data
                :class:`plotly.graph_objects.layout.template.Da
                ta` instance or dict with compatible properties
            layout
                :class:`plotly.graph_objects.Layout` instance
                or dict with compatible properties
""",
            ),
            **kwargs,
        )
