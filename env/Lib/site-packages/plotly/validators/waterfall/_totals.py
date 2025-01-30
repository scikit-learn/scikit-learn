import _plotly_utils.basevalidators


class TotalsValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="totals", parent_name="waterfall", **kwargs):
        super(TotalsValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Totals"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            marker
                :class:`plotly.graph_objects.waterfall.totals.M
                arker` instance or dict with compatible
                properties
""",
            ),
            **kwargs,
        )
