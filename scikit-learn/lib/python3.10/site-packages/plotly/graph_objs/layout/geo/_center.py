#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Center(_BaseLayoutHierarchyType):
    _parent_path_str = "layout.geo"
    _path_str = "layout.geo.center"
    _valid_props = {"lat", "lon"}

    @property
    def lat(self):
        """
        Sets the latitude of the map's center. For all projection
        types, the map's latitude center lies at the middle of the
        latitude range by default.

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

    @property
    def lon(self):
        """
        Sets the longitude of the map's center. By default, the map's
        longitude center lies at the middle of the longitude range for
        scoped projection and above `projection.rotation.lon`
        otherwise.

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

    @property
    def _prop_descriptions(self):
        return """\
        lat
            Sets the latitude of the map's center. For all
            projection types, the map's latitude center lies at the
            middle of the latitude range by default.
        lon
            Sets the longitude of the map's center. By default, the
            map's longitude center lies at the middle of the
            longitude range for scoped projection and above
            `projection.rotation.lon` otherwise.
        """

    def __init__(self, arg=None, lat=None, lon=None, **kwargs):
        """
        Construct a new Center object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.layout.geo.Center`
        lat
            Sets the latitude of the map's center. For all
            projection types, the map's latitude center lies at the
            middle of the latitude range by default.
        lon
            Sets the longitude of the map's center. By default, the
            map's longitude center lies at the middle of the
            longitude range for scoped projection and above
            `projection.rotation.lon` otherwise.

        Returns
        -------
        Center
        """
        super().__init__("center")
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
The first argument to the plotly.graph_objs.layout.geo.Center
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.geo.Center`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("lat", arg, lat)
        self._set_property("lon", arg, lon)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
