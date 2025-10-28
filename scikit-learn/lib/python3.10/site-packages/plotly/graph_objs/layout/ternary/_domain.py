#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Domain(_BaseLayoutHierarchyType):
    _parent_path_str = "layout.ternary"
    _path_str = "layout.ternary.domain"
    _valid_props = {"column", "row", "x", "y"}

    @property
    def column(self):
        """
        If there is a layout grid, use the domain for this column in
        the grid for this ternary subplot .

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

    @property
    def row(self):
        """
        If there is a layout grid, use the domain for this row in the
        grid for this ternary subplot .

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

    @property
    def x(self):
        """
            Sets the horizontal domain of this ternary subplot (in plot
            fraction).

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

    @property
    def y(self):
        """
            Sets the vertical domain of this ternary subplot (in plot
            fraction).

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

    @property
    def _prop_descriptions(self):
        return """\
        column
            If there is a layout grid, use the domain for this
            column in the grid for this ternary subplot .
        row
            If there is a layout grid, use the domain for this row
            in the grid for this ternary subplot .
        x
            Sets the horizontal domain of this ternary subplot (in
            plot fraction).
        y
            Sets the vertical domain of this ternary subplot (in
            plot fraction).
        """

    def __init__(self, arg=None, column=None, row=None, x=None, y=None, **kwargs):
        """
        Construct a new Domain object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.layout.ternary.Domain`
        column
            If there is a layout grid, use the domain for this
            column in the grid for this ternary subplot .
        row
            If there is a layout grid, use the domain for this row
            in the grid for this ternary subplot .
        x
            Sets the horizontal domain of this ternary subplot (in
            plot fraction).
        y
            Sets the vertical domain of this ternary subplot (in
            plot fraction).

        Returns
        -------
        Domain
        """
        super().__init__("domain")
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
The first argument to the plotly.graph_objs.layout.ternary.Domain
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.ternary.Domain`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("column", arg, column)
        self._set_property("row", arg, row)
        self._set_property("x", arg, x)
        self._set_property("y", arg, y)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
