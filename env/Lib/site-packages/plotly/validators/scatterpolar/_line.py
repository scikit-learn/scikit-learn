import _plotly_utils.basevalidators


class LineValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="line", parent_name="scatterpolar", **kwargs):
        super(LineValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Line"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            backoff
                Sets the line back off from the end point of
                the nth line segment (in px). This option is
                useful e.g. to avoid overlap with arrowhead
                markers. With "auto" the lines would trim
                before markers if `marker.angleref` is set to
                "previous".
            backoffsrc
                Sets the source reference on Chart Studio Cloud
                for `backoff`.
            color
                Sets the line color.
            dash
                Sets the dash style of lines. Set to a dash
                type string ("solid", "dot", "dash",
                "longdash", "dashdot", or "longdashdot") or a
                dash length list in px (eg "5px,10px,2px,2px").
            shape
                Determines the line shape. With "spline" the
                lines are drawn using spline interpolation. The
                other available values correspond to step-wise
                line shapes.
            smoothing
                Has an effect only if `shape` is set to
                "spline" Sets the amount of smoothing. 0
                corresponds to no smoothing (equivalent to a
                "linear" shape).
            width
                Sets the line width (in px).
""",
            ),
            **kwargs,
        )
