from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Domain(_BaseLayoutHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "layout.geo"
    _path_str = "layout.geo.domain"
    _valid_props = {"column", "row", "x", "y"}

    # column
    # ------
    @property
    def column(self):
        """
        If there is a layout grid, use the domain for this column in
        the grid for this geo subplot . Note that geo subplots are
        constrained by domain. In general, when `projection.scale` is
        set to 1. a map will fit either its x or y domain, but not
        both.

        The 'column' property is a integer and may be specified as:
          - An int (or float that will be cast to an int)
            in the interval [0, 9223372036854775807]

        Returns
        -------
        int
        """
        return self["column"]

    @column.setter
    def column(self, val):
        self["column"] = val

    # row
    # ---
    @property
    def row(self):
        """
        If there is a layout grid, use the domain for this row in the
        grid for this geo subplot . Note that geo subplots are
        constrained by domain. In general, when `projection.scale` is
        set to 1. a map will fit either its x or y domain, but not
        both.

        The 'row' property is a integer and may be specified as:
          - An int (or float that will be cast to an int)
            in the interval [0, 9223372036854775807]

        Returns
        -------
        int
        """
        return self["row"]

    @row.setter
    def row(self, val):
        self["row"] = val

    # x
    # -
    @property
    def x(self):
        """
            Sets the horizontal domain of this geo subplot (in plot
            fraction). Note that geo subplots are constrained by domain. In
            general, when `projection.scale` is set to 1. a map will fit
            either its x or y domain, but not both.

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
            Sets the vertical domain of this geo subplot (in plot
            fraction). Note that geo subplots are constrained by domain. In
            general, when `projection.scale` is set to 1. a map will fit
            either its x or y domain, but not both.

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
        column
            If there is a layout grid, use the domain for this
            column in the grid for this geo subplot . Note that geo
            subplots are constrained by domain. In general, when
            `projection.scale` is set to 1. a map will fit either
            its x or y domain, but not both.
        row
            If there is a layout grid, use the domain for this row
            in the grid for this geo subplot . Note that geo
            subplots are constrained by domain. In general, when
            `projection.scale` is set to 1. a map will fit either
            its x or y domain, but not both.
        x
            Sets the horizontal domain of this geo subplot (in plot
            fraction). Note that geo subplots are constrained by
            domain. In general, when `projection.scale` is set to
            1. a map will fit either its x or y domain, but not
            both.
        y
            Sets the vertical domain of this geo subplot (in plot
            fraction). Note that geo subplots are constrained by
            domain. In general, when `projection.scale` is set to
            1. a map will fit either its x or y domain, but not
            both.
        """

    def __init__(self, arg=None, column=None, row=None, x=None, y=None, **kwargs):
        """
        Construct a new Domain object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.layout.geo.Domain`
        column
            If there is a layout grid, use the domain for this
            column in the grid for this geo subplot . Note that geo
            subplots are constrained by domain. In general, when
            `projection.scale` is set to 1. a map will fit either
            its x or y domain, but not both.
        row
            If there is a layout grid, use the domain for this row
            in the grid for this geo subplot . Note that geo
            subplots are constrained by domain. In general, when
            `projection.scale` is set to 1. a map will fit either
            its x or y domain, but not both.
        x
            Sets the horizontal domain of this geo subplot (in plot
            fraction). Note that geo subplots are constrained by
            domain. In general, when `projection.scale` is set to
            1. a map will fit either its x or y domain, but not
            both.
        y
            Sets the vertical domain of this geo subplot (in plot
            fraction). Note that geo subplots are constrained by
            domain. In general, when `projection.scale` is set to
            1. a map will fit either its x or y domain, but not
            both.

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
The first argument to the plotly.graph_objs.layout.geo.Domain
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.geo.Domain`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("column", None)
        _v = column if column is not None else _v
        if _v is not None:
            self["column"] = _v
        _v = arg.pop("row", None)
        _v = row if row is not None else _v
        if _v is not None:
            self["row"] = _v
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
