from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class YAxis(_BaseLayoutHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "layout.xaxis.rangeslider"
    _path_str = "layout.xaxis.rangeslider.yaxis"
    _valid_props = {"range", "rangemode"}

    # range
    # -----
    @property
    def range(self):
        """
            Sets the range of this axis for the rangeslider.

            The 'range' property is an info array that may be specified as:

            * a list or tuple of 2 elements where:
        (0) The 'range[0]' property accepts values of any type
        (1) The 'range[1]' property accepts values of any type

            Returns
            -------
            list
        """
        return self["range"]

    @range.setter
    def range(self, val):
        self["range"] = val

    # rangemode
    # ---------
    @property
    def rangemode(self):
        """
        Determines whether or not the range of this axis in the
        rangeslider use the same value than in the main plot when
        zooming in/out. If "auto", the autorange will be used. If
        "fixed", the `range` is used. If "match", the current range of
        the corresponding y-axis on the main subplot is used.

        The 'rangemode' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['auto', 'fixed', 'match']

        Returns
        -------
        Any
        """
        return self["rangemode"]

    @rangemode.setter
    def rangemode(self, val):
        self["rangemode"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        range
            Sets the range of this axis for the rangeslider.
        rangemode
            Determines whether or not the range of this axis in the
            rangeslider use the same value than in the main plot
            when zooming in/out. If "auto", the autorange will be
            used. If "fixed", the `range` is used. If "match", the
            current range of the corresponding y-axis on the main
            subplot is used.
        """

    def __init__(self, arg=None, range=None, rangemode=None, **kwargs):
        """
        Construct a new YAxis object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of :class:`plotly.graph_objs.layout.xaxis.r
            angeslider.YAxis`
        range
            Sets the range of this axis for the rangeslider.
        rangemode
            Determines whether or not the range of this axis in the
            rangeslider use the same value than in the main plot
            when zooming in/out. If "auto", the autorange will be
            used. If "fixed", the `range` is used. If "match", the
            current range of the corresponding y-axis on the main
            subplot is used.

        Returns
        -------
        YAxis
        """
        super(YAxis, self).__init__("yaxis")

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
The first argument to the plotly.graph_objs.layout.xaxis.rangeslider.YAxis
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.xaxis.rangeslider.YAxis`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("range", None)
        _v = range if range is not None else _v
        if _v is not None:
            self["range"] = _v
        _v = arg.pop("rangemode", None)
        _v = rangemode if rangemode is not None else _v
        if _v is not None:
            self["rangemode"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
