from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Rangebreak(_BaseLayoutHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "layout.xaxis"
    _path_str = "layout.xaxis.rangebreak"
    _valid_props = {
        "bounds",
        "dvalue",
        "enabled",
        "name",
        "pattern",
        "templateitemname",
        "values",
    }

    # bounds
    # ------
    @property
    def bounds(self):
        """
            Sets the lower and upper bounds of this axis rangebreak. Can be
            used with `pattern`.

            The 'bounds' property is an info array that may be specified as:

            * a list or tuple of 2 elements where:
        (0) The 'bounds[0]' property accepts values of any type
        (1) The 'bounds[1]' property accepts values of any type

            Returns
            -------
            list
        """
        return self["bounds"]

    @bounds.setter
    def bounds(self, val):
        self["bounds"] = val

    # dvalue
    # ------
    @property
    def dvalue(self):
        """
        Sets the size of each `values` item. The default is one day in
        milliseconds.

        The 'dvalue' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["dvalue"]

    @dvalue.setter
    def dvalue(self, val):
        self["dvalue"] = val

    # enabled
    # -------
    @property
    def enabled(self):
        """
        Determines whether this axis rangebreak is enabled or disabled.
        Please note that `rangebreaks` only work for "date" axis type.

        The 'enabled' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["enabled"]

    @enabled.setter
    def enabled(self, val):
        self["enabled"] = val

    # name
    # ----
    @property
    def name(self):
        """
        When used in a template, named items are created in the output
        figure in addition to any items the figure already has in this
        array. You can modify these items in the output figure by
        making your own item with `templateitemname` matching this
        `name` alongside your modifications (including `visible: false`
        or `enabled: false` to hide it). Has no effect outside of a
        template.

        The 'name' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["name"]

    @name.setter
    def name(self, val):
        self["name"] = val

    # pattern
    # -------
    @property
    def pattern(self):
        """
        Determines a pattern on the time line that generates breaks. If
        *day of week* - days of the week in English e.g. 'Sunday' or
        `sun` (matching is case-insensitive and considers only the
        first three characters), as well as Sunday-based integers
        between 0 and 6. If "hour" - hour (24-hour clock) as decimal
        numbers between 0 and 24. for more info. Examples: - { pattern:
        'day of week', bounds: [6, 1] }  or simply { bounds: ['sat',
        'mon'] }   breaks from Saturday to Monday (i.e. skips the
        weekends). - { pattern: 'hour', bounds: [17, 8] }   breaks from
        5pm to 8am (i.e. skips non-work hours).

        The 'pattern' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['day of week', 'hour', '']

        Returns
        -------
        Any
        """
        return self["pattern"]

    @pattern.setter
    def pattern(self, val):
        self["pattern"] = val

    # templateitemname
    # ----------------
    @property
    def templateitemname(self):
        """
        Used to refer to a named item in this array in the template.
        Named items from the template will be created even without a
        matching item in the input figure, but you can modify one by
        making an item with `templateitemname` matching its `name`,
        alongside your modifications (including `visible: false` or
        `enabled: false` to hide it). If there is no template or no
        matching item, this item will be hidden unless you explicitly
        show it with `visible: true`.

        The 'templateitemname' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["templateitemname"]

    @templateitemname.setter
    def templateitemname(self, val):
        self["templateitemname"] = val

    # values
    # ------
    @property
    def values(self):
        """
        Sets the coordinate values corresponding to the rangebreaks. An
        alternative to `bounds`. Use `dvalue` to set the size of the
        values along the axis.

        The 'values' property is an info array that may be specified as:
        * a list of elements where:
          The 'values[i]' property accepts values of any type

        Returns
        -------
        list
        """
        return self["values"]

    @values.setter
    def values(self, val):
        self["values"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        bounds
            Sets the lower and upper bounds of this axis
            rangebreak. Can be used with `pattern`.
        dvalue
            Sets the size of each `values` item. The default is one
            day in milliseconds.
        enabled
            Determines whether this axis rangebreak is enabled or
            disabled. Please note that `rangebreaks` only work for
            "date" axis type.
        name
            When used in a template, named items are created in the
            output figure in addition to any items the figure
            already has in this array. You can modify these items
            in the output figure by making your own item with
            `templateitemname` matching this `name` alongside your
            modifications (including `visible: false` or `enabled:
            false` to hide it). Has no effect outside of a
            template.
        pattern
            Determines a pattern on the time line that generates
            breaks. If *day of week* - days of the week in English
            e.g. 'Sunday' or `sun` (matching is case-insensitive
            and considers only the first three characters), as well
            as Sunday-based integers between 0 and 6. If "hour" -
            hour (24-hour clock) as decimal numbers between 0 and
            24. for more info. Examples: - { pattern: 'day of
            week', bounds: [6, 1] }  or simply { bounds: ['sat',
            'mon'] }   breaks from Saturday to Monday (i.e. skips
            the weekends). - { pattern: 'hour', bounds: [17, 8] }
            breaks from 5pm to 8am (i.e. skips non-work hours).
        templateitemname
            Used to refer to a named item in this array in the
            template. Named items from the template will be created
            even without a matching item in the input figure, but
            you can modify one by making an item with
            `templateitemname` matching its `name`, alongside your
            modifications (including `visible: false` or `enabled:
            false` to hide it). If there is no template or no
            matching item, this item will be hidden unless you
            explicitly show it with `visible: true`.
        values
            Sets the coordinate values corresponding to the
            rangebreaks. An alternative to `bounds`. Use `dvalue`
            to set the size of the values along the axis.
        """

    def __init__(
        self,
        arg=None,
        bounds=None,
        dvalue=None,
        enabled=None,
        name=None,
        pattern=None,
        templateitemname=None,
        values=None,
        **kwargs,
    ):
        """
        Construct a new Rangebreak object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.layout.xaxis.Rangebreak`
        bounds
            Sets the lower and upper bounds of this axis
            rangebreak. Can be used with `pattern`.
        dvalue
            Sets the size of each `values` item. The default is one
            day in milliseconds.
        enabled
            Determines whether this axis rangebreak is enabled or
            disabled. Please note that `rangebreaks` only work for
            "date" axis type.
        name
            When used in a template, named items are created in the
            output figure in addition to any items the figure
            already has in this array. You can modify these items
            in the output figure by making your own item with
            `templateitemname` matching this `name` alongside your
            modifications (including `visible: false` or `enabled:
            false` to hide it). Has no effect outside of a
            template.
        pattern
            Determines a pattern on the time line that generates
            breaks. If *day of week* - days of the week in English
            e.g. 'Sunday' or `sun` (matching is case-insensitive
            and considers only the first three characters), as well
            as Sunday-based integers between 0 and 6. If "hour" -
            hour (24-hour clock) as decimal numbers between 0 and
            24. for more info. Examples: - { pattern: 'day of
            week', bounds: [6, 1] }  or simply { bounds: ['sat',
            'mon'] }   breaks from Saturday to Monday (i.e. skips
            the weekends). - { pattern: 'hour', bounds: [17, 8] }
            breaks from 5pm to 8am (i.e. skips non-work hours).
        templateitemname
            Used to refer to a named item in this array in the
            template. Named items from the template will be created
            even without a matching item in the input figure, but
            you can modify one by making an item with
            `templateitemname` matching its `name`, alongside your
            modifications (including `visible: false` or `enabled:
            false` to hide it). If there is no template or no
            matching item, this item will be hidden unless you
            explicitly show it with `visible: true`.
        values
            Sets the coordinate values corresponding to the
            rangebreaks. An alternative to `bounds`. Use `dvalue`
            to set the size of the values along the axis.

        Returns
        -------
        Rangebreak
        """
        super(Rangebreak, self).__init__("rangebreaks")

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
The first argument to the plotly.graph_objs.layout.xaxis.Rangebreak
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.xaxis.Rangebreak`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("bounds", None)
        _v = bounds if bounds is not None else _v
        if _v is not None:
            self["bounds"] = _v
        _v = arg.pop("dvalue", None)
        _v = dvalue if dvalue is not None else _v
        if _v is not None:
            self["dvalue"] = _v
        _v = arg.pop("enabled", None)
        _v = enabled if enabled is not None else _v
        if _v is not None:
            self["enabled"] = _v
        _v = arg.pop("name", None)
        _v = name if name is not None else _v
        if _v is not None:
            self["name"] = _v
        _v = arg.pop("pattern", None)
        _v = pattern if pattern is not None else _v
        if _v is not None:
            self["pattern"] = _v
        _v = arg.pop("templateitemname", None)
        _v = templateitemname if templateitemname is not None else _v
        if _v is not None:
            self["templateitemname"] = _v
        _v = arg.pop("values", None)
        _v = values if values is not None else _v
        if _v is not None:
            self["values"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
