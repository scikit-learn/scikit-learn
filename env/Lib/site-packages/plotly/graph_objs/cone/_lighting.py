from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Lighting(_BaseTraceHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "cone"
    _path_str = "cone.lighting"
    _valid_props = {
        "ambient",
        "diffuse",
        "facenormalsepsilon",
        "fresnel",
        "roughness",
        "specular",
        "vertexnormalsepsilon",
    }

    # ambient
    # -------
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

    # diffuse
    # -------
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

    # facenormalsepsilon
    # ------------------
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

    # fresnel
    # -------
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

    # roughness
    # ---------
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

    # specular
    # --------
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

    # vertexnormalsepsilon
    # --------------------
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

    # Self properties description
    # ---------------------------
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
            an instance of :class:`plotly.graph_objs.cone.Lighting`
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
        super(Lighting, self).__init__("lighting")

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
The first argument to the plotly.graph_objs.cone.Lighting
constructor must be a dict or
an instance of :class:`plotly.graph_objs.cone.Lighting`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("ambient", None)
        _v = ambient if ambient is not None else _v
        if _v is not None:
            self["ambient"] = _v
        _v = arg.pop("diffuse", None)
        _v = diffuse if diffuse is not None else _v
        if _v is not None:
            self["diffuse"] = _v
        _v = arg.pop("facenormalsepsilon", None)
        _v = facenormalsepsilon if facenormalsepsilon is not None else _v
        if _v is not None:
            self["facenormalsepsilon"] = _v
        _v = arg.pop("fresnel", None)
        _v = fresnel if fresnel is not None else _v
        if _v is not None:
            self["fresnel"] = _v
        _v = arg.pop("roughness", None)
        _v = roughness if roughness is not None else _v
        if _v is not None:
            self["roughness"] = _v
        _v = arg.pop("specular", None)
        _v = specular if specular is not None else _v
        if _v is not None:
            self["specular"] = _v
        _v = arg.pop("vertexnormalsepsilon", None)
        _v = vertexnormalsepsilon if vertexnormalsepsilon is not None else _v
        if _v is not None:
            self["vertexnormalsepsilon"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
