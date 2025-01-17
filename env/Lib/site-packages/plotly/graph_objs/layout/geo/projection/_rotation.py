from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Rotation(_BaseLayoutHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "layout.geo.projection"
    _path_str = "layout.geo.projection.rotation"
    _valid_props = {"lat", "lon", "roll"}

    # lat
    # ---
    @property
    def lat(self):
        """
        Rotates the map along meridians (in degrees North).

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
        Rotates the map along parallels (in degrees East). Defaults to
        the center of the `lonaxis.range` values.

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

    # roll
    # ----
    @property
    def roll(self):
        """
        Roll the map (in degrees) For example, a roll of 180 makes the
        map appear upside down.

        The 'roll' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["roll"]

    @roll.setter
    def roll(self, val):
        self["roll"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        lat
            Rotates the map along meridians (in degrees North).
        lon
            Rotates the map along parallels (in degrees East).
            Defaults to the center of the `lonaxis.range` values.
        roll
            Roll the map (in degrees) For example, a roll of 180
            makes the map appear upside down.
        """

    def __init__(self, arg=None, lat=None, lon=None, roll=None, **kwargs):
        """
        Construct a new Rotation object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of :class:`plotly.graph_objs.layout.geo.pro
            jection.Rotation`
        lat
            Rotates the map along meridians (in degrees North).
        lon
            Rotates the map along parallels (in degrees East).
            Defaults to the center of the `lonaxis.range` values.
        roll
            Roll the map (in degrees) For example, a roll of 180
            makes the map appear upside down.

        Returns
        -------
        Rotation
        """
        super(Rotation, self).__init__("rotation")

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
The first argument to the plotly.graph_objs.layout.geo.projection.Rotation
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.geo.projection.Rotation`"""
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
        _v = arg.pop("roll", None)
        _v = roll if roll is not None else _v
        if _v is not None:
            self["roll"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
