#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Lighting(_BaseTraceHierarchyType):
    _parent_path_str = "isosurface"
    _path_str = "isosurface.lighting"
    _valid_props = {
        "ambient",
        "diffuse",
        "facenormalsepsilon",
        "fresnel",
        "roughness",
        "specular",
        "vertexnormalsepsilon",
    }

    @property
    def ambient(self):
        """
        Ambient light increases overall color visibility but can wash
        out the image.

        The 'ambient' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["ambient"]

    @ambient.setter
    def ambient(self, val):
        self["ambient"] = val

    @property
    def diffuse(self):
        """
        Represents the extent that incident rays are reflected in a
        range of angles.

        The 'diffuse' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["diffuse"]

    @diffuse.setter
    def diffuse(self, val):
        self["diffuse"] = val

    @property
    def facenormalsepsilon(self):
        """
        Epsilon for face normals calculation avoids math issues arising
        from degenerate geometry.

        The 'facenormalsepsilon' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["facenormalsepsilon"]

    @facenormalsepsilon.setter
    def facenormalsepsilon(self, val):
        self["facenormalsepsilon"] = val

    @property
    def fresnel(self):
        """
        Represents the reflectance as a dependency of the viewing
        angle; e.g. paper is reflective when viewing it from the edge
        of the paper (almost 90 degrees), causing shine.

        The 'fresnel' property is a number and may be specified as:
          - An int or float in the interval [0, 5]

        Returns
        -------
        int|float
        """
        return self["fresnel"]

    @fresnel.setter
    def fresnel(self, val):
        self["fresnel"] = val

    @property
    def roughness(self):
        """
        Alters specular reflection; the rougher the surface, the wider
        and less contrasty the shine.

        The 'roughness' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["roughness"]

    @roughness.setter
    def roughness(self, val):
        self["roughness"] = val

    @property
    def specular(self):
        """
        Represents the level that incident rays are reflected in a
        single direction, causing shine.

        The 'specular' property is a number and may be specified as:
          - An int or float in the interval [0, 2]

        Returns
        -------
        int|float
        """
        return self["specular"]

    @specular.setter
    def specular(self, val):
        self["specular"] = val

    @property
    def vertexnormalsepsilon(self):
        """
        Epsilon for vertex normals calculation avoids math issues
        arising from degenerate geometry.

        The 'vertexnormalsepsilon' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["vertexnormalsepsilon"]

    @vertexnormalsepsilon.setter
    def vertexnormalsepsilon(self, val):
        self["vertexnormalsepsilon"] = val

    @property
    def _prop_descriptions(self):
        return """\
        ambient
            Ambient light increases overall color visibility but
            can wash out the image.
        diffuse
            Represents the extent that incident rays are reflected
            in a range of angles.
        facenormalsepsilon
            Epsilon for face normals calculation avoids math issues
            arising from degenerate geometry.
        fresnel
            Represents the reflectance as a dependency of the
            viewing angle; e.g. paper is reflective when viewing it
            from the edge of the paper (almost 90 degrees), causing
            shine.
        roughness
            Alters specular reflection; the rougher the surface,
            the wider and less contrasty the shine.
        specular
            Represents the level that incident rays are reflected
            in a single direction, causing shine.
        vertexnormalsepsilon
            Epsilon for vertex normals calculation avoids math
            issues arising from degenerate geometry.
        """

    def __init__(
        self,
        arg=None,
        ambient=None,
        diffuse=None,
        facenormalsepsilon=None,
        fresnel=None,
        roughness=None,
        specular=None,
        vertexnormalsepsilon=None,
        **kwargs,
    ):
        """
        Construct a new Lighting object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.isosurface.Lighting`
        ambient
            Ambient light increases overall color visibility but
            can wash out the image.
        diffuse
            Represents the extent that incident rays are reflected
            in a range of angles.
        facenormalsepsilon
            Epsilon for face normals calculation avoids math issues
            arising from degenerate geometry.
        fresnel
            Represents the reflectance as a dependency of the
            viewing angle; e.g. paper is reflective when viewing it
            from the edge of the paper (almost 90 degrees), causing
            shine.
        roughness
            Alters specular reflection; the rougher the surface,
            the wider and less contrasty the shine.
        specular
            Represents the level that incident rays are reflected
            in a single direction, causing shine.
        vertexnormalsepsilon
            Epsilon for vertex normals calculation avoids math
            issues arising from degenerate geometry.

        Returns
        -------
        Lighting
        """
        super().__init__("lighting")
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
The first argument to the plotly.graph_objs.isosurface.Lighting
constructor must be a dict or
an instance of :class:`plotly.graph_objs.isosurface.Lighting`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("ambient", arg, ambient)
        self._set_property("diffuse", arg, diffuse)
        self._set_property("facenormalsepsilon", arg, facenormalsepsilon)
        self._set_property("fresnel", arg, fresnel)
        self._set_property("roughness", arg, roughness)
        self._set_property("specular", arg, specular)
        self._set_property("vertexnormalsepsilon", arg, vertexnormalsepsilon)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
