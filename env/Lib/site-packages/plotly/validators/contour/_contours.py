import _plotly_utils.basevalidators


class ContoursValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="contours", parent_name="contour", **kwargs):
        super(ContoursValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Contours"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            coloring
                Determines the coloring method showing the
                contour values. If "fill", coloring is done
                evenly between each contour level If "heatmap",
                a heatmap gradient coloring is applied between
                each contour level. If "lines", coloring is
                done on the contour lines. If "none", no
                coloring is applied on this trace.
            end
                Sets the end contour level value. Must be more
                than `contours.start`
            labelfont
                Sets the font used for labeling the contour
                levels. The default color comes from the lines,
                if shown. The default family and size come from
                `layout.font`.
            labelformat
                Sets the contour label formatting rule using d3
                formatting mini-languages which are very
                similar to those in Python. For numbers, see: h
                ttps://github.com/d3/d3-format/tree/v1.4.5#d3-
                format.
            operation
                Sets the constraint operation. "=" keeps
                regions equal to `value` "<" and "<=" keep
                regions less than `value` ">" and ">=" keep
                regions greater than `value` "[]", "()", "[)",
                and "(]" keep regions inside `value[0]` to
                `value[1]` "][", ")(", "](", ")[" keep regions
                outside `value[0]` to value[1]` Open vs. closed
                intervals make no difference to constraint
                display, but all versions are allowed for
                consistency with filter transforms.
            showlabels
                Determines whether to label the contour lines
                with their values.
            showlines
                Determines whether or not the contour lines are
                drawn. Has an effect only if
                `contours.coloring` is set to "fill".
            size
                Sets the step between each contour level. Must
                be positive.
            start
                Sets the starting contour level value. Must be
                less than `contours.end`
            type
                If `levels`, the data is represented as a
                contour plot with multiple levels displayed. If
                `constraint`, the data is represented as
                constraints with the invalid region shaded as
                specified by the `operation` and `value`
                parameters.
            value
                Sets the value or values of the constraint
                boundary. When `operation` is set to one of the
                comparison values (=,<,>=,>,<=) "value" is
                expected to be a number. When `operation` is
                set to one of the interval values
                ([],(),[),(],][,)(,](,)[) "value" is expected
                to be an array of two numbers where the first
                is the lower bound and the second is the upper
                bound.
""",
            ),
            **kwargs,
        )
