#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Projection(_BaseLayoutHierarchyType):
    _parent_path_str = "layout.geo"
    _path_str = "layout.geo.projection"
    _valid_props = {"distance", "parallels", "rotation", "scale", "tilt", "type"}

    @property
    def distance(self):
        """
        For satellite projection type only. Sets the distance from the
        center of the sphere to the point of view as a proportion of
        the sphere’s radius.

        The 'distance' property is a number and may be specified as:
          - An int or float in the interval [1.001, inf]

        Returns
        -------
        int|float
        """
        return self["distance"]

    @distance.setter
    def distance(self, val):
        self["distance"] = val

    @property
    def parallels(self):
        """
            For conic projection types only. Sets the parallels (tangent,
            secant) where the cone intersects the sphere.

            The 'parallels' property is an info array that may be specified as:

            * a list or tuple of 2 elements where:
        (0) The 'parallels[0]' property is a number and may be specified as:
              - An int or float
        (1) The 'parallels[1]' property is a number and may be specified as:
              - An int or float

            Returns
            -------
            list
        """
        return self["parallels"]

    @parallels.setter
    def parallels(self, val):
        self["parallels"] = val

    @property
    def rotation(self):
        """
        The 'rotation' property is an instance of Rotation
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.geo.projection.Rotation`
          - A dict of string/value properties that will be passed
            to the Rotation constructor

        Returns
        -------
        plotly.graph_objs.layout.geo.projection.Rotation
        """
        return self["rotation"]

    @rotation.setter
    def rotation(self, val):
        self["rotation"] = val

    @property
    def scale(self):
        """
        Zooms in or out on the map view. A scale of 1 corresponds to
        the largest zoom level that fits the map's lon and lat ranges.

        The 'scale' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["scale"]

    @scale.setter
    def scale(self, val):
        self["scale"] = val

    @property
    def tilt(self):
        """
        For satellite projection type only. Sets the tilt angle of
        perspective projection.

        The 'tilt' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["tilt"]

    @tilt.setter
    def tilt(self, val):
        self["tilt"] = val

    @property
    def type(self):
        """
        Sets the projection type.

        The 'type' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['airy', 'aitoff', 'albers', 'albers usa', 'august',
                'azimuthal equal area', 'azimuthal equidistant', 'baker',
                'bertin1953', 'boggs', 'bonne', 'bottomley', 'bromley',
                'collignon', 'conic conformal', 'conic equal area', 'conic
                equidistant', 'craig', 'craster', 'cylindrical equal
                area', 'cylindrical stereographic', 'eckert1', 'eckert2',
                'eckert3', 'eckert4', 'eckert5', 'eckert6', 'eisenlohr',
                'equal earth', 'equirectangular', 'fahey', 'foucaut',
                'foucaut sinusoidal', 'ginzburg4', 'ginzburg5',
                'ginzburg6', 'ginzburg8', 'ginzburg9', 'gnomonic',
                'gringorten', 'gringorten quincuncial', 'guyou', 'hammer',
                'hill', 'homolosine', 'hufnagel', 'hyperelliptical',
                'kavrayskiy7', 'lagrange', 'larrivee', 'laskowski',
                'loximuthal', 'mercator', 'miller', 'mollweide', 'mt flat
                polar parabolic', 'mt flat polar quartic', 'mt flat polar
                sinusoidal', 'natural earth', 'natural earth1', 'natural
                earth2', 'nell hammer', 'nicolosi', 'orthographic',
                'patterson', 'peirce quincuncial', 'polyconic',
                'rectangular polyconic', 'robinson', 'satellite', 'sinu
                mollweide', 'sinusoidal', 'stereographic', 'times',
                'transverse mercator', 'van der grinten', 'van der
                grinten2', 'van der grinten3', 'van der grinten4',
                'wagner4', 'wagner6', 'wiechel', 'winkel tripel',
                'winkel3']

        Returns
        -------
        Any
        """
        return self["type"]

    @type.setter
    def type(self, val):
        self["type"] = val

    @property
    def _prop_descriptions(self):
        return """\
        distance
            For satellite projection type only. Sets the distance
            from the center of the sphere to the point of view as a
            proportion of the sphere’s radius.
        parallels
            For conic projection types only. Sets the parallels
            (tangent, secant) where the cone intersects the sphere.
        rotation
            :class:`plotly.graph_objects.layout.geo.projection.Rota
            tion` instance or dict with compatible properties
        scale
            Zooms in or out on the map view. A scale of 1
            corresponds to the largest zoom level that fits the
            map's lon and lat ranges.
        tilt
            For satellite projection type only. Sets the tilt angle
            of perspective projection.
        type
            Sets the projection type.
        """

    def __init__(
        self,
        arg=None,
        distance=None,
        parallels=None,
        rotation=None,
        scale=None,
        tilt=None,
        type=None,
        **kwargs,
    ):
        """
        Construct a new Projection object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.layout.geo.Projection`
        distance
            For satellite projection type only. Sets the distance
            from the center of the sphere to the point of view as a
            proportion of the sphere’s radius.
        parallels
            For conic projection types only. Sets the parallels
            (tangent, secant) where the cone intersects the sphere.
        rotation
            :class:`plotly.graph_objects.layout.geo.projection.Rota
            tion` instance or dict with compatible properties
        scale
            Zooms in or out on the map view. A scale of 1
            corresponds to the largest zoom level that fits the
            map's lon and lat ranges.
        tilt
            For satellite projection type only. Sets the tilt angle
            of perspective projection.
        type
            Sets the projection type.

        Returns
        -------
        Projection
        """
        super().__init__("projection")
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
The first argument to the plotly.graph_objs.layout.geo.Projection
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.geo.Projection`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("distance", arg, distance)
        self._set_property("parallels", arg, parallels)
        self._set_property("rotation", arg, rotation)
        self._set_property("scale", arg, scale)
        self._set_property("tilt", arg, tilt)
        self._set_property("type", arg, type)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
