from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Totals(_BaseTraceHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "waterfall"
    _path_str = "waterfall.totals"
    _valid_props = {"marker"}

    # marker
    # ------
    @property
    def marker(self):
        """
        The 'marker' property is an instance of Marker
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.waterfall.totals.Marker`
          - A dict of string/value properties that will be passed
            to the Marker constructor

            Supported dict properties:

                color
                    Sets the marker color of all intermediate sums
                    and total values.
                line
                    :class:`plotly.graph_objects.waterfall.totals.m
                    arker.Line` instance or dict with compatible
                    properties

        Returns
        -------
        plotly.graph_objs.waterfall.totals.Marker
        """
        return self["marker"]

    @marker.setter
    def marker(self, val):
        self["marker"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        marker
            :class:`plotly.graph_objects.waterfall.totals.Marker`
            instance or dict with compatible properties
        """

    def __init__(self, arg=None, marker=None, **kwargs):
        """
        Construct a new Totals object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.waterfall.Totals`
        marker
            :class:`plotly.graph_objects.waterfall.totals.Marker`
            instance or dict with compatible properties

        Returns
        -------
        Totals
        """
        super(Totals, self).__init__("totals")

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
The first argument to the plotly.graph_objs.waterfall.Totals
constructor must be a dict or
an instance of :class:`plotly.graph_objs.waterfall.Totals`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("marker", None)
        _v = marker if marker is not None else _v
        if _v is not None:
            self["marker"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
