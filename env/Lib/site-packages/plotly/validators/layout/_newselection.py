import _plotly_utils.basevalidators


class NewselectionValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="newselection", parent_name="layout", **kwargs):
        super(NewselectionValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Newselection"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            line
                :class:`plotly.graph_objects.layout.newselectio
                n.Line` instance or dict with compatible
                properties
            mode
                Describes how a new selection is created. If
                `immediate`, a new selection is created after
                first mouse up. If `gradual`, a new selection
                is not created after first mouse. By adding to
                and subtracting from the initial selection,
                this option allows declaring extra outlines of
                the selection.
""",
            ),
            **kwargs,
        )
