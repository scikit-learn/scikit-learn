import _plotly_utils.basevalidators


class TernaryValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="ternary", parent_name="layout", **kwargs):
        super(TernaryValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Ternary"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            aaxis
                :class:`plotly.graph_objects.layout.ternary.Aax
                is` instance or dict with compatible properties
            baxis
                :class:`plotly.graph_objects.layout.ternary.Bax
                is` instance or dict with compatible properties
            bgcolor
                Set the background color of the subplot
            caxis
                :class:`plotly.graph_objects.layout.ternary.Cax
                is` instance or dict with compatible properties
            domain
                :class:`plotly.graph_objects.layout.ternary.Dom
                ain` instance or dict with compatible
                properties
            sum
                The number each triplet should sum to, and the
                maximum range of each axis
            uirevision
                Controls persistence of user-driven changes in
                axis `min` and `title`, if not overridden in
                the individual axes. Defaults to
                `layout.uirevision`.
""",
            ),
            **kwargs,
        )
