#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Bounds(_BaseLayoutHierarchyType):
    _parent_path_str = "layout.mapbox"
    _path_str = "layout.mapbox.bounds"
    _valid_props = {"east", "north", "south", "west"}

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
        super().__init__("bounds")
        if "_parent" in kwargs:
            self._parent = kwargs["_parent"]
            return

        if arg is None:
            arg = {}
        elif isinstance(arg, self.__class__):
            arg = arg.to_plotly_json()
        elif isinstance(arg, dict):
            arg = _copy.copy(arg)
        else:
            raise ValueError("""\
The first argument to the plotly.graph_objs.layout.mapbox.Bounds
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.mapbox.Bounds`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("east", arg, east)
        self._set_property("north", arg, north)
        self._set_property("south", arg, south)
        self._set_property("west", arg, west)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
