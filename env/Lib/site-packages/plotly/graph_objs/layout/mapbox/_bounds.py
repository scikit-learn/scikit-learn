from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Bounds(_BaseLayoutHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "layout.mapbox"
    _path_str = "layout.mapbox.bounds"
    _valid_props = {"east", "north", "south", "west"}

    # east
    # ----
    @property
    def east(self):
        """
        Sets the maximum longitude of the map (in degrees East) if
        `west`, `south` and `north` are declared.

        The 'east' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["east"]

    @east.setter
    def east(self, val):
        self["east"] = val

    # north
    # -----
    @property
    def north(self):
        """
        Sets the maximum latitude of the map (in degrees North) if
        `east`, `west` and `south` are declared.

        The 'north' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["north"]

    @north.setter
    def north(self, val):
        self["north"] = val

    # south
    # -----
    @property
    def south(self):
        """
        Sets the minimum latitude of the map (in degrees North) if
        `east`, `west` and `north` are declared.

        The 'south' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["south"]

    @south.setter
    def south(self, val):
        self["south"] = val

    # west
    # ----
    @property
    def west(self):
        """
        Sets the minimum longitude of the map (in degrees East) if
        `east`, `south` and `north` are declared.

        The 'west' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["west"]

    @west.setter
    def west(self, val):
        self["west"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        east
            Sets the maximum longitude of the map (in degrees East)
            if `west`, `south` and `north` are declared.
        north
            Sets the maximum latitude of the map (in degrees North)
            if `east`, `west` and `south` are declared.
        south
            Sets the minimum latitude of the map (in degrees North)
            if `east`, `west` and `north` are declared.
        west
            Sets the minimum longitude of the map (in degrees East)
            if `east`, `south` and `north` are declared.
        """

    def __init__(
        self, arg=None, east=None, north=None, south=None, west=None, **kwargs
    ):
        """
        Construct a new Bounds object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.layout.mapbox.Bounds`
        east
            Sets the maximum longitude of the map (in degrees East)
            if `west`, `south` and `north` are declared.
        north
            Sets the maximum latitude of the map (in degrees North)
            if `east`, `west` and `south` are declared.
        south
            Sets the minimum latitude of the map (in degrees North)
            if `east`, `west` and `north` are declared.
        west
            Sets the minimum longitude of the map (in degrees East)
            if `east`, `south` and `north` are declared.

        Returns
        -------
        Bounds
        """
        super(Bounds, self).__init__("bounds")

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
The first argument to the plotly.graph_objs.layout.mapbox.Bounds
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.mapbox.Bounds`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("east", None)
        _v = east if east is not None else _v
        if _v is not None:
            self["east"] = _v
        _v = arg.pop("north", None)
        _v = north if north is not None else _v
        if _v is not None:
            self["north"] = _v
        _v = arg.pop("south", None)
        _v = south if south is not None else _v
        if _v is not None:
            self["south"] = _v
        _v = arg.pop("west", None)
        _v = west if west is not None else _v
        if _v is not None:
            self["west"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
