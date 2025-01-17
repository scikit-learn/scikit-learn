import _plotly_utils.basevalidators


class DomainValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="domain", parent_name="layout.geo", **kwargs):
        super(DomainValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Domain"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            column
                If there is a layout grid, use the domain for
                this column in the grid for this geo subplot .
                Note that geo subplots are constrained by
                domain. In general, when `projection.scale` is
                set to 1. a map will fit either its x or y
                domain, but not both.
            row
                If there is a layout grid, use the domain for
                this row in the grid for this geo subplot .
                Note that geo subplots are constrained by
                domain. In general, when `projection.scale` is
                set to 1. a map will fit either its x or y
                domain, but not both.
            x
                Sets the horizontal domain of this geo subplot
                (in plot fraction). Note that geo subplots are
                constrained by domain. In general, when
                `projection.scale` is set to 1. a map will fit
                either its x or y domain, but not both.
            y
                Sets the vertical domain of this geo subplot
                (in plot fraction). Note that geo subplots are
                constrained by domain. In general, when
                `projection.scale` is set to 1. a map will fit
                either its x or y domain, but not both.
""",
            ),
            **kwargs,
        )
