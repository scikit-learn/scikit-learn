from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Projection(_BaseTraceHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "scatter3d"
    _path_str = "scatter3d.projection"
    _valid_props = {"x", "y", "z"}

    # x
    # -
    @property
    def x(self):
        """
        The 'x' property is an instance of X
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.scatter3d.projection.X`
          - A dict of string/value properties that will be passed
            to the X constructor

            Supported dict properties:

                opacity
                    Sets the projection color.
                scale
                    Sets the scale factor determining the size of
                    the projection marker points.
                show
                    Sets whether or not projections are shown along
                    the x axis.

        Returns
        -------
        plotly.graph_objs.scatter3d.projection.X
        """
        return self["x"]

    @x.setter
    def x(self, val):
        self["x"] = val

    # y
    # -
    @property
    def y(self):
        """
        The 'y' property is an instance of Y
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.scatter3d.projection.Y`
          - A dict of string/value properties that will be passed
            to the Y constructor

            Supported dict properties:

                opacity
                    Sets the projection color.
                scale
                    Sets the scale factor determining the size of
                    the projection marker points.
                show
                    Sets whether or not projections are shown along
                    the y axis.

        Returns
        -------
        plotly.graph_objs.scatter3d.projection.Y
        """
        return self["y"]

    @y.setter
    def y(self, val):
        self["y"] = val

    # z
    # -
    @property
    def z(self):
        """
        The 'z' property is an instance of Z
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.scatter3d.projection.Z`
          - A dict of string/value properties that will be passed
            to the Z constructor

            Supported dict properties:

                opacity
                    Sets the projection color.
                scale
                    Sets the scale factor determining the size of
                    the projection marker points.
                show
                    Sets whether or not projections are shown along
                    the z axis.

        Returns
        -------
        plotly.graph_objs.scatter3d.projection.Z
        """
        return self["z"]

    @z.setter
    def z(self, val):
        self["z"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        x
            :class:`plotly.graph_objects.scatter3d.projection.X`
            instance or dict with compatible properties
        y
            :class:`plotly.graph_objects.scatter3d.projection.Y`
            instance or dict with compatible properties
        z
            :class:`plotly.graph_objects.scatter3d.projection.Z`
            instance or dict with compatible properties
        """

    def __init__(self, arg=None, x=None, y=None, z=None, **kwargs):
        """
        Construct a new Projection object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.scatter3d.Projection`
        x
            :class:`plotly.graph_objects.scatter3d.projection.X`
            instance or dict with compatible properties
        y
            :class:`plotly.graph_objects.scatter3d.projection.Y`
            instance or dict with compatible properties
        z
            :class:`plotly.graph_objects.scatter3d.projection.Z`
            instance or dict with compatible properties

        Returns
        -------
        Projection
        """
        super(Projection, self).__init__("projection")

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
The first argument to the plotly.graph_objs.scatter3d.Projection
constructor must be a dict or
an instance of :class:`plotly.graph_objs.scatter3d.Projection`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("x", None)
        _v = x if x is not None else _v
        if _v is not None:
            self["x"] = _v
        _v = arg.pop("y", None)
        _v = y if y is not None else _v
        if _v is not None:
            self["y"] = _v
        _v = arg.pop("z", None)
        _v = z if z is not None else _v
        if _v is not None:
            self["z"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
