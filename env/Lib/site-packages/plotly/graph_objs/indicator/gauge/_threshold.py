from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Threshold(_BaseTraceHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "indicator.gauge"
    _path_str = "indicator.gauge.threshold"
    _valid_props = {"line", "thickness", "value"}

    # line
    # ----
    @property
    def line(self):
        """
        The 'line' property is an instance of Line
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.indicator.gauge.threshold.Line`
          - A dict of string/value properties that will be passed
            to the Line constructor

            Supported dict properties:

                color
                    Sets the color of the threshold line.
                width
                    Sets the width (in px) of the threshold line.

        Returns
        -------
        plotly.graph_objs.indicator.gauge.threshold.Line
        """
        return self["line"]

    @line.setter
    def line(self, val):
        self["line"] = val

    # thickness
    # ---------
    @property
    def thickness(self):
        """
        Sets the thickness of the threshold line as a fraction of the
        thickness of the gauge.

        The 'thickness' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["thickness"]

    @thickness.setter
    def thickness(self, val):
        self["thickness"] = val

    # value
    # -----
    @property
    def value(self):
        """
        Sets a treshold value drawn as a line.

        The 'value' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["value"]

    @value.setter
    def value(self, val):
        self["value"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        line
            :class:`plotly.graph_objects.indicator.gauge.threshold.
            Line` instance or dict with compatible properties
        thickness
            Sets the thickness of the threshold line as a fraction
            of the thickness of the gauge.
        value
            Sets a treshold value drawn as a line.
        """

    def __init__(self, arg=None, line=None, thickness=None, value=None, **kwargs):
        """
        Construct a new Threshold object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.indicator.gauge.Threshold`
        line
            :class:`plotly.graph_objects.indicator.gauge.threshold.
            Line` instance or dict with compatible properties
        thickness
            Sets the thickness of the threshold line as a fraction
            of the thickness of the gauge.
        value
            Sets a treshold value drawn as a line.

        Returns
        -------
        Threshold
        """
        super(Threshold, self).__init__("threshold")

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
The first argument to the plotly.graph_objs.indicator.gauge.Threshold
constructor must be a dict or
an instance of :class:`plotly.graph_objs.indicator.gauge.Threshold`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("line", None)
        _v = line if line is not None else _v
        if _v is not None:
            self["line"] = _v
        _v = arg.pop("thickness", None)
        _v = thickness if thickness is not None else _v
        if _v is not None:
            self["thickness"] = _v
        _v = arg.pop("value", None)
        _v = value if value is not None else _v
        if _v is not None:
            self["value"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
