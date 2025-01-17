import _plotly_utils.basevalidators


class YbinsValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="ybins", parent_name="histogram2d", **kwargs):
        super(YbinsValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "YBins"),
            data_docs=kwargs.pop(
                "data_docs",
                """
            end
                Sets the end value for the y axis bins. The
                last bin may not end exactly at this value, we
                increment the bin edge by `size` from `start`
                until we reach or exceed `end`. Defaults to the
                maximum data value. Like `start`, for dates use
                a date string, and for category data `end` is
                based on the category serial numbers.
            size
                Sets the size of each y axis bin. Default
                behavior: If `nbinsy` is 0 or omitted, we
                choose a nice round bin size such that the
                number of bins is about the same as the typical
                number of samples in each bin. If `nbinsy` is
                provided, we choose a nice round bin size
                giving no more than that many bins. For date
                data, use milliseconds or "M<n>" for months, as
                in `axis.dtick`. For category data, the number
                of categories to bin together (always defaults
                to 1).
            start
                Sets the starting value for the y axis bins.
                Defaults to the minimum data value, shifted
                down if necessary to make nice round values and
                to remove ambiguous bin edges. For example, if
                most of the data is integers we shift the bin
                edges 0.5 down, so a `size` of 5 would have a
                default `start` of -0.5, so it is clear that
                0-4 are in the first bin, 5-9 in the second,
                but continuous data gets a start of 0 and bins
                [0,5), [5,10) etc. Dates behave similarly, and
                `start` should be a date string. For category
                data, `start` is based on the category serial
                numbers, and defaults to -0.5.
""",
            ),
            **kwargs,
        )
