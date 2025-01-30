import _plotly_utils.basevalidators


class CellsValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="cells", parent_name="table", **kwargs):
        super(CellsValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Cells"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            align
                Sets the horizontal alignment of the `text`
                within the box. Has an effect only if `text`
                spans two or more lines (i.e. `text` contains
                one or more <br> HTML tags) or if an explicit
                width is set to override the text width.
            alignsrc
                Sets the source reference on Chart Studio Cloud
                for `align`.
            fill
                :class:`plotly.graph_objects.table.cells.Fill`
                instance or dict with compatible properties
            font
                :class:`plotly.graph_objects.table.cells.Font`
                instance or dict with compatible properties
            format
                Sets the cell value formatting rule using d3
                formatting mini-languages which are very
                similar to those in Python. For numbers, see: h
                ttps://github.com/d3/d3-format/tree/v1.4.5#d3-
                format.
            formatsrc
                Sets the source reference on Chart Studio Cloud
                for `format`.
            height
                The height of cells.
            line
                :class:`plotly.graph_objects.table.cells.Line`
                instance or dict with compatible properties
            prefix
                Prefix for cell values.
            prefixsrc
                Sets the source reference on Chart Studio Cloud
                for `prefix`.
            suffix
                Suffix for cell values.
            suffixsrc
                Sets the source reference on Chart Studio Cloud
                for `suffix`.
            values
                Cell values. `values[m][n]` represents the
                value of the `n`th point in column `m`,
                therefore the `values[m]` vector length for all
                columns must be the same (longer vectors will
                be truncated). Each value must be a finite
                number or a string.
            valuessrc
                Sets the source reference on Chart Studio Cloud
                for `values`.
""",
            ),
            **kwargs,
        )
