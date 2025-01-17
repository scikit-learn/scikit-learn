import _plotly_utils.basevalidators


class DomainValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="domain", parent_name="layout.grid", **kwargs):
        super(DomainValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Domain"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            x
                Sets the horizontal domain of this grid subplot
                (in plot fraction). The first and last cells
                end exactly at the domain edges, with no grout
                around the edges.
            y
                Sets the vertical domain of this grid subplot
                (in plot fraction). The first and last cells
                end exactly at the domain edges, with no grout
                around the edges.
""",
            ),
            **kwargs,
        )
