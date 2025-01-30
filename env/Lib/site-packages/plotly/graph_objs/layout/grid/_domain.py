from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Domain(_BaseLayoutHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "layout.grid"
    _path_str = "layout.grid.domain"
    _valid_props = {"x", "y"}

    # x
    # -
    @property
    def x(self):
        """
            Sets the horizontal domain of this grid subplot (in plot
            fraction). The first and last cells end exactly at the domain
            edges, with no grout around the edges.

            The 'x' property is an info array that may be specified as:

            * a list or tuple of 2 elements where:
        (0) The 'x[0]' property is a number and may be specified as:
              - An int or float in the interval [0, 1]
        (1) The 'x[1]' property is a number and may be specified as:
              - An int or float in the interval [0, 1]

            Returns
            -------
            list
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
            Sets the vertical domain of this grid subplot (in plot
            fraction). The first and last cells end exactly at the domain
            edges, with no grout around the edges.

            The 'y' property is an info array that may be specified as:

            * a list or tuple of 2 elements where:
        (0) The 'y[0]' property is a number and may be specified as:
              - An int or float in the interval [0, 1]
        (1) The 'y[1]' property is a number and may be specified as:
              - An int or float in the interval [0, 1]

            Returns
            -------
            list
        """
        return self["y"]

    @y.setter
    def y(self, val):
        self["y"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        x
            Sets the horizontal domain of this grid subplot (in
            plot fraction). The first and last cells end exactly at
            the domain edges, with no grout around the edges.
        y
            Sets the vertical domain of this grid subplot (in plot
            fraction). The first and last cells end exactly at the
            domain edges, with no grout around the edges.
        """

    def __init__(self, arg=None, x=None, y=None, **kwargs):
        """
        Construct a new Domain object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.layout.grid.Domain`
        x
            Sets the horizontal domain of this grid subplot (in
            plot fraction). The first and last cells end exactly at
            the domain edges, with no grout around the edges.
        y
            Sets the vertical domain of this grid subplot (in plot
            fraction). The first and last cells end exactly at the
            domain edges, with no grout around the edges.

        Returns
        -------
        Domain
        """
        super(Domain, self).__init__("domain")

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
The first argument to the plotly.graph_objs.layout.grid.Domain
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.grid.Domain`"""
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

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
