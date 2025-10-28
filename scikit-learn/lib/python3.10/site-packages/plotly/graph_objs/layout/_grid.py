#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Grid(_BaseLayoutHierarchyType):
    _parent_path_str = "layout"
    _path_str = "layout.grid"
    _valid_props = {
        "columns",
        "domain",
        "pattern",
        "roworder",
        "rows",
        "subplots",
        "xaxes",
        "xgap",
        "xside",
        "yaxes",
        "ygap",
        "yside",
    }

    @property
    def columns(self):
        """
        The number of columns in the grid. If you provide a 2D
        `subplots` array, the length of its longest row is used as the
        default. If you give an `xaxes` array, its length is used as
        the default. But it's also possible to have a different length,
        if you want to leave a row at the end for non-cartesian
        subplots.

        The 'columns' property is a integer and may be specified as:
          - An int (or float that will be cast to an int)
            in the interval [1, 9223372036854775807]

        Returns
        -------
        int
        """
        return self["columns"]

    @columns.setter
    def columns(self, val):
        self["columns"] = val

    @property
    def domain(self):
        """
        The 'domain' property is an instance of Domain
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.grid.Domain`
          - A dict of string/value properties that will be passed
            to the Domain constructor

        Returns
        -------
        plotly.graph_objs.layout.grid.Domain
        """
        return self["domain"]

    @domain.setter
    def domain(self, val):
        self["domain"] = val

    @property
    def pattern(self):
        """
        If no `subplots`, `xaxes`, or `yaxes` are given but we do have
        `rows` and `columns`, we can generate defaults using
        consecutive axis IDs, in two ways: "coupled" gives one x axis
        per column and one y axis per row. "independent" uses a new xy
        pair for each cell, left-to-right across each row then
        iterating rows according to `roworder`.

        The 'pattern' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['independent', 'coupled']

        Returns
        -------
        Any
        """
        return self["pattern"]

    @pattern.setter
    def pattern(self, val):
        self["pattern"] = val

    @property
    def roworder(self):
        """
        Is the first row the top or the bottom? Note that columns are
        always enumerated from left to right.

        The 'roworder' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['top to bottom', 'bottom to top']

        Returns
        -------
        Any
        """
        return self["roworder"]

    @roworder.setter
    def roworder(self, val):
        self["roworder"] = val

    @property
    def rows(self):
        """
        The number of rows in the grid. If you provide a 2D `subplots`
        array or a `yaxes` array, its length is used as the default.
        But it's also possible to have a different length, if you want
        to leave a row at the end for non-cartesian subplots.

        The 'rows' property is a integer and may be specified as:
          - An int (or float that will be cast to an int)
            in the interval [1, 9223372036854775807]

        Returns
        -------
        int
        """
        return self["rows"]

    @rows.setter
    def rows(self, val):
        self["rows"] = val

    @property
    def subplots(self):
        """
        Used for freeform grids, where some axes may be shared across
        subplots but others are not. Each entry should be a cartesian
        subplot id, like "xy" or "x3y2", or "" to leave that cell
        empty. You may reuse x axes within the same column, and y axes
        within the same row. Non-cartesian subplots and traces that
        support `domain` can place themselves in this grid separately
        using the `gridcell` attribute.

        The 'subplots' property is an info array that may be specified as:
        * a 2D list where:
          The 'subplots[i][j]' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['']
          - A string that matches one of the following regular expressions:
                ['^x([2-9]|[1-9][0-9]+)?y([2-9]|[1-9][0-9]+)?$']

        Returns
        -------
        list
        """
        return self["subplots"]

    @subplots.setter
    def subplots(self, val):
        self["subplots"] = val

    @property
    def xaxes(self):
        """
        Used with `yaxes` when the x and y axes are shared across
        columns and rows. Each entry should be an x axis id like "x",
        "x2", etc., or "" to not put an x axis in that column. Entries
        other than "" must be unique. Ignored if `subplots` is present.
        If missing but `yaxes` is present, will generate consecutive
        IDs.

        The 'xaxes' property is an info array that may be specified as:
        * a list of elements where:
          The 'xaxes[i]' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['']
          - A string that matches one of the following regular expressions:
                ['^x([2-9]|[1-9][0-9]+)?( domain)?$']

        Returns
        -------
        list
        """
        return self["xaxes"]

    @xaxes.setter
    def xaxes(self, val):
        self["xaxes"] = val

    @property
    def xgap(self):
        """
        Horizontal space between grid cells, expressed as a fraction of
        the total width available to one cell. Defaults to 0.1 for
        coupled-axes grids and 0.2 for independent grids.

        The 'xgap' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["xgap"]

    @xgap.setter
    def xgap(self, val):
        self["xgap"] = val

    @property
    def xside(self):
        """
        Sets where the x axis labels and titles go. "bottom" means the
        very bottom of the grid. "bottom plot" is the lowest plot that
        each x axis is used in. "top" and "top plot" are similar.

        The 'xside' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['bottom', 'bottom plot', 'top plot', 'top']

        Returns
        -------
        Any
        """
        return self["xside"]

    @xside.setter
    def xside(self, val):
        self["xside"] = val

    @property
    def yaxes(self):
        """
        Used with `yaxes` when the x and y axes are shared across
        columns and rows. Each entry should be an y axis id like "y",
        "y2", etc., or "" to not put a y axis in that row. Entries
        other than "" must be unique. Ignored if `subplots` is present.
        If missing but `xaxes` is present, will generate consecutive
        IDs.

        The 'yaxes' property is an info array that may be specified as:
        * a list of elements where:
          The 'yaxes[i]' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['']
          - A string that matches one of the following regular expressions:
                ['^y([2-9]|[1-9][0-9]+)?( domain)?$']

        Returns
        -------
        list
        """
        return self["yaxes"]

    @yaxes.setter
    def yaxes(self, val):
        self["yaxes"] = val

    @property
    def ygap(self):
        """
        Vertical space between grid cells, expressed as a fraction of
        the total height available to one cell. Defaults to 0.1 for
        coupled-axes grids and 0.3 for independent grids.

        The 'ygap' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["ygap"]

    @ygap.setter
    def ygap(self, val):
        self["ygap"] = val

    @property
    def yside(self):
        """
        Sets where the y axis labels and titles go. "left" means the
        very left edge of the grid. *left plot* is the leftmost plot
        that each y axis is used in. "right" and *right plot* are
        similar.

        The 'yside' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['left', 'left plot', 'right plot', 'right']

        Returns
        -------
        Any
        """
        return self["yside"]

    @yside.setter
    def yside(self, val):
        self["yside"] = val

    @property
    def _prop_descriptions(self):
        return """\
        columns
            The number of columns in the grid. If you provide a 2D
            `subplots` array, the length of its longest row is used
            as the default. If you give an `xaxes` array, its
            length is used as the default. But it's also possible
            to have a different length, if you want to leave a row
            at the end for non-cartesian subplots.
        domain
            :class:`plotly.graph_objects.layout.grid.Domain`
            instance or dict with compatible properties
        pattern
            If no `subplots`, `xaxes`, or `yaxes` are given but we
            do have `rows` and `columns`, we can generate defaults
            using consecutive axis IDs, in two ways: "coupled"
            gives one x axis per column and one y axis per row.
            "independent" uses a new xy pair for each cell, left-
            to-right across each row then iterating rows according
            to `roworder`.
        roworder
            Is the first row the top or the bottom? Note that
            columns are always enumerated from left to right.
        rows
            The number of rows in the grid. If you provide a 2D
            `subplots` array or a `yaxes` array, its length is used
            as the default. But it's also possible to have a
            different length, if you want to leave a row at the end
            for non-cartesian subplots.
        subplots
            Used for freeform grids, where some axes may be shared
            across subplots but others are not. Each entry should
            be a cartesian subplot id, like "xy" or "x3y2", or ""
            to leave that cell empty. You may reuse x axes within
            the same column, and y axes within the same row. Non-
            cartesian subplots and traces that support `domain` can
            place themselves in this grid separately using the
            `gridcell` attribute.
        xaxes
            Used with `yaxes` when the x and y axes are shared
            across columns and rows. Each entry should be an x axis
            id like "x", "x2", etc., or "" to not put an x axis in
            that column. Entries other than "" must be unique.
            Ignored if `subplots` is present. If missing but
            `yaxes` is present, will generate consecutive IDs.
        xgap
            Horizontal space between grid cells, expressed as a
            fraction of the total width available to one cell.
            Defaults to 0.1 for coupled-axes grids and 0.2 for
            independent grids.
        xside
            Sets where the x axis labels and titles go. "bottom"
            means the very bottom of the grid. "bottom plot" is the
            lowest plot that each x axis is used in. "top" and "top
            plot" are similar.
        yaxes
            Used with `yaxes` when the x and y axes are shared
            across columns and rows. Each entry should be an y axis
            id like "y", "y2", etc., or "" to not put a y axis in
            that row. Entries other than "" must be unique. Ignored
            if `subplots` is present. If missing but `xaxes` is
            present, will generate consecutive IDs.
        ygap
            Vertical space between grid cells, expressed as a
            fraction of the total height available to one cell.
            Defaults to 0.1 for coupled-axes grids and 0.3 for
            independent grids.
        yside
            Sets where the y axis labels and titles go. "left"
            means the very left edge of the grid. *left plot* is
            the leftmost plot that each y axis is used in. "right"
            and *right plot* are similar.
        """

    def __init__(
        self,
        arg=None,
        columns=None,
        domain=None,
        pattern=None,
        roworder=None,
        rows=None,
        subplots=None,
        xaxes=None,
        xgap=None,
        xside=None,
        yaxes=None,
        ygap=None,
        yside=None,
        **kwargs,
    ):
        """
        Construct a new Grid object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of :class:`plotly.graph_objs.layout.Grid`
        columns
            The number of columns in the grid. If you provide a 2D
            `subplots` array, the length of its longest row is used
            as the default. If you give an `xaxes` array, its
            length is used as the default. But it's also possible
            to have a different length, if you want to leave a row
            at the end for non-cartesian subplots.
        domain
            :class:`plotly.graph_objects.layout.grid.Domain`
            instance or dict with compatible properties
        pattern
            If no `subplots`, `xaxes`, or `yaxes` are given but we
            do have `rows` and `columns`, we can generate defaults
            using consecutive axis IDs, in two ways: "coupled"
            gives one x axis per column and one y axis per row.
            "independent" uses a new xy pair for each cell, left-
            to-right across each row then iterating rows according
            to `roworder`.
        roworder
            Is the first row the top or the bottom? Note that
            columns are always enumerated from left to right.
        rows
            The number of rows in the grid. If you provide a 2D
            `subplots` array or a `yaxes` array, its length is used
            as the default. But it's also possible to have a
            different length, if you want to leave a row at the end
            for non-cartesian subplots.
        subplots
            Used for freeform grids, where some axes may be shared
            across subplots but others are not. Each entry should
            be a cartesian subplot id, like "xy" or "x3y2", or ""
            to leave that cell empty. You may reuse x axes within
            the same column, and y axes within the same row. Non-
            cartesian subplots and traces that support `domain` can
            place themselves in this grid separately using the
            `gridcell` attribute.
        xaxes
            Used with `yaxes` when the x and y axes are shared
            across columns and rows. Each entry should be an x axis
            id like "x", "x2", etc., or "" to not put an x axis in
            that column. Entries other than "" must be unique.
            Ignored if `subplots` is present. If missing but
            `yaxes` is present, will generate consecutive IDs.
        xgap
            Horizontal space between grid cells, expressed as a
            fraction of the total width available to one cell.
            Defaults to 0.1 for coupled-axes grids and 0.2 for
            independent grids.
        xside
            Sets where the x axis labels and titles go. "bottom"
            means the very bottom of the grid. "bottom plot" is the
            lowest plot that each x axis is used in. "top" and "top
            plot" are similar.
        yaxes
            Used with `yaxes` when the x and y axes are shared
            across columns and rows. Each entry should be an y axis
            id like "y", "y2", etc., or "" to not put a y axis in
            that row. Entries other than "" must be unique. Ignored
            if `subplots` is present. If missing but `xaxes` is
            present, will generate consecutive IDs.
        ygap
            Vertical space between grid cells, expressed as a
            fraction of the total height available to one cell.
            Defaults to 0.1 for coupled-axes grids and 0.3 for
            independent grids.
        yside
            Sets where the y axis labels and titles go. "left"
            means the very left edge of the grid. *left plot* is
            the leftmost plot that each y axis is used in. "right"
            and *right plot* are similar.

        Returns
        -------
        Grid
        """
        super().__init__("grid")
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
The first argument to the plotly.graph_objs.layout.Grid
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.Grid`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("columns", arg, columns)
        self._set_property("domain", arg, domain)
        self._set_property("pattern", arg, pattern)
        self._set_property("roworder", arg, roworder)
        self._set_property("rows", arg, rows)
        self._set_property("subplots", arg, subplots)
        self._set_property("xaxes", arg, xaxes)
        self._set_property("xgap", arg, xgap)
        self._set_property("xside", arg, xside)
        self._set_property("yaxes", arg, yaxes)
        self._set_property("ygap", arg, ygap)
        self._set_property("yside", arg, yside)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
