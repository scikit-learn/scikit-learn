from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Center(_BaseLayoutHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "layout.mapbox"
    _path_str = "layout.mapbox.center"
    _valid_props = {"lat", "lon"}

    # lat
    # ---
    @property
    def lat(self):
        """
        Sets the latitude of the center of the map (in degrees North).

        The 'lat' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["lat"]

    @lat.setter
    def lat(self, val):
        self["lat"] = val

    # lon
    # ---
    @property
    def lon(self):
        """
        Sets the longitude of the center of the map (in degrees East).

        The 'lon' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["lon"]

    @lon.setter
    def lon(self, val):
        self["lon"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        lat
            Sets the latitude of the center of the map (in degrees
            North).
        lon
            Sets the longitude of the center of the map (in degrees
            East).
        """

    def __init__(self, arg=None, lat=None, lon=None, **kwargs):
        """
        Construct a new Center object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.layout.mapbox.Center`
        lat
            Sets the latitude of the center of the map (in degrees
            North).
        lon
            Sets the longitude of the center of the map (in degrees
            East).

        Returns
        -------
        Center
        """
        super(Center, self).__init__("center")

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
The first argument to the plotly.graph_objs.layout.mapbox.Center
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.mapbox.Center`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("lat", None)
        _v = lat if lat is not None else _v
        if _v is not None:
            self["lat"] = _v
        _v = arg.pop("lon", None)
        _v = lon if lon is not None else _v
        if _v is not None:
            self["lon"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
