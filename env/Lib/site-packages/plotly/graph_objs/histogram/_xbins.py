from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class XBins(_BaseTraceHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "histogram"
    _path_str = "histogram.xbins"
    _valid_props = {"end", "size", "start"}

    # end
    # ---
    @property
    def end(self):
        """
        Sets the end value for the x axis bins. The last bin may not
        end exactly at this value, we increment the bin edge by `size`
        from `start` until we reach or exceed `end`. Defaults to the
        maximum data value. Like `start`, for dates use a date string,
        and for category data `end` is based on the category serial
        numbers.

        The 'end' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["end"]

    @end.setter
    def end(self, val):
        self["end"] = val

    # size
    # ----
    @property
    def size(self):
        """
        Sets the size of each x axis bin. Default behavior: If `nbinsx`
        is 0 or omitted, we choose a nice round bin size such that the
        number of bins is about the same as the typical number of
        samples in each bin. If `nbinsx` is provided, we choose a nice
        round bin size giving no more than that many bins. For date
        data, use milliseconds or "M<n>" for months, as in
        `axis.dtick`. For category data, the number of categories to
        bin together (always defaults to 1). If multiple non-overlaying
        histograms share a subplot, the first explicit `size` is used
        and all others discarded. If no `size` is provided,the sample
        data from all traces is combined to determine `size` as
        described above.

        The 'size' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["size"]

    @size.setter
    def size(self, val):
        self["size"] = val

    # start
    # -----
    @property
    def start(self):
        """
        Sets the starting value for the x axis bins. Defaults to the
        minimum data value, shifted down if necessary to make nice
        round values and to remove ambiguous bin edges. For example, if
        most of the data is integers we shift the bin edges 0.5 down,
        so a `size` of 5 would have a default `start` of -0.5, so it is
        clear that 0-4 are in the first bin, 5-9 in the second, but
        continuous data gets a start of 0 and bins [0,5), [5,10) etc.
        Dates behave similarly, and `start` should be a date string.
        For category data, `start` is based on the category serial
        numbers, and defaults to -0.5. If multiple non-overlaying
        histograms share a subplot, the first explicit `start` is used
        exactly and all others are shifted down (if necessary) to
        differ from that one by an integer number of bins.

        The 'start' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["start"]

    @start.setter
    def start(self, val):
        self["start"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        end
            Sets the end value for the x axis bins. The last bin
            may not end exactly at this value, we increment the bin
            edge by `size` from `start` until we reach or exceed
            `end`. Defaults to the maximum data value. Like
            `start`, for dates use a date string, and for category
            data `end` is based on the category serial numbers.
        size
            Sets the size of each x axis bin. Default behavior: If
            `nbinsx` is 0 or omitted, we choose a nice round bin
            size such that the number of bins is about the same as
            the typical number of samples in each bin. If `nbinsx`
            is provided, we choose a nice round bin size giving no
            more than that many bins. For date data, use
            milliseconds or "M<n>" for months, as in `axis.dtick`.
            For category data, the number of categories to bin
            together (always defaults to 1). If multiple non-
            overlaying histograms share a subplot, the first
            explicit `size` is used and all others discarded. If no
            `size` is provided,the sample data from all traces is
            combined to determine `size` as described above.
        start
            Sets the starting value for the x axis bins. Defaults
            to the minimum data value, shifted down if necessary to
            make nice round values and to remove ambiguous bin
            edges. For example, if most of the data is integers we
            shift the bin edges 0.5 down, so a `size` of 5 would
            have a default `start` of -0.5, so it is clear that 0-4
            are in the first bin, 5-9 in the second, but continuous
            data gets a start of 0 and bins [0,5), [5,10) etc.
            Dates behave similarly, and `start` should be a date
            string. For category data, `start` is based on the
            category serial numbers, and defaults to -0.5. If
            multiple non-overlaying histograms share a subplot, the
            first explicit `start` is used exactly and all others
            are shifted down (if necessary) to differ from that one
            by an integer number of bins.
        """

    def __init__(self, arg=None, end=None, size=None, start=None, **kwargs):
        """
        Construct a new XBins object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.histogram.XBins`
        end
            Sets the end value for the x axis bins. The last bin
            may not end exactly at this value, we increment the bin
            edge by `size` from `start` until we reach or exceed
            `end`. Defaults to the maximum data value. Like
            `start`, for dates use a date string, and for category
            data `end` is based on the category serial numbers.
        size
            Sets the size of each x axis bin. Default behavior: If
            `nbinsx` is 0 or omitted, we choose a nice round bin
            size such that the number of bins is about the same as
            the typical number of samples in each bin. If `nbinsx`
            is provided, we choose a nice round bin size giving no
            more than that many bins. For date data, use
            milliseconds or "M<n>" for months, as in `axis.dtick`.
            For category data, the number of categories to bin
            together (always defaults to 1). If multiple non-
            overlaying histograms share a subplot, the first
            explicit `size` is used and all others discarded. If no
            `size` is provided,the sample data from all traces is
            combined to determine `size` as described above.
        start
            Sets the starting value for the x axis bins. Defaults
            to the minimum data value, shifted down if necessary to
            make nice round values and to remove ambiguous bin
            edges. For example, if most of the data is integers we
            shift the bin edges 0.5 down, so a `size` of 5 would
            have a default `start` of -0.5, so it is clear that 0-4
            are in the first bin, 5-9 in the second, but continuous
            data gets a start of 0 and bins [0,5), [5,10) etc.
            Dates behave similarly, and `start` should be a date
            string. For category data, `start` is based on the
            category serial numbers, and defaults to -0.5. If
            multiple non-overlaying histograms share a subplot, the
            first explicit `start` is used exactly and all others
            are shifted down (if necessary) to differ from that one
            by an integer number of bins.

        Returns
        -------
        XBins
        """
        super(XBins, self).__init__("xbins")

        if "_parent" in kwargs:
            self._parent = kwargs["_parent"]
            return

        # Validate arg
        # ------------
        if arg is None:
            arg = {}
        elif isinstance(arg, self.__class__):
            arg = arg.to_plotly_json()
        elif isinstance(arg, dict):
            arg = _copy.copy(arg)
        else:
            raise ValueError(
                """\
The first argument to the plotly.graph_objs.histogram.XBins
constructor must be a dict or
an instance of :class:`plotly.graph_objs.histogram.XBins`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("end", None)
        _v = end if end is not None else _v
        if _v is not None:
            self["end"] = _v
        _v = arg.pop("size", None)
        _v = size if size is not None else _v
        if _v is not None:
            self["size"] = _v
        _v = arg.pop("start", None)
        _v = start if start is not None else _v
        if _v is not None:
            self["start"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
