from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Button(_BaseLayoutHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "layout.xaxis.rangeselector"
    _path_str = "layout.xaxis.rangeselector.button"
    _valid_props = {
        "count",
        "label",
        "name",
        "step",
        "stepmode",
        "templateitemname",
        "visible",
    }

    # count
    # -----
    @property
    def count(self):
        """
        Sets the number of steps to take to update the range. Use with
        `step` to specify the update interval.

        The 'count' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["count"]

    @count.setter
    def count(self, val):
        self["count"] = val

    # label
    # -----
    @property
    def label(self):
        """
        Sets the text label to appear on the button.

        The 'label' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["label"]

    @label.setter
    def label(self, val):
        self["label"] = val

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

    # step
    # ----
    @property
    def step(self):
        """
        The unit of measurement that the `count` value will set the
        range by.

        The 'step' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['month', 'year', 'day', 'hour', 'minute', 'second',
                'all']

        Returns
        -------
        Any
        """
        return self["step"]

    @step.setter
    def step(self, val):
        self["step"] = val

    # stepmode
    # --------
    @property
    def stepmode(self):
        """
        Sets the range update mode. If "backward", the range update
        shifts the start of range back "count" times "step"
        milliseconds. If "todate", the range update shifts the start of
        range back to the first timestamp from "count" times "step"
        milliseconds back. For example, with `step` set to "year" and
        `count` set to 1 the range update shifts the start of the range
        back to January 01 of the current year. Month and year "todate"
        are currently available only for the built-in (Gregorian)
        calendar.

        The 'stepmode' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['backward', 'todate']

        Returns
        -------
        Any
        """
        return self["stepmode"]

    @stepmode.setter
    def stepmode(self, val):
        self["stepmode"] = val

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

    # visible
    # -------
    @property
    def visible(self):
        """
        Determines whether or not this button is visible.

        The 'visible' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["visible"]

    @visible.setter
    def visible(self, val):
        self["visible"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        count
            Sets the number of steps to take to update the range.
            Use with `step` to specify the update interval.
        label
            Sets the text label to appear on the button.
        name
            When used in a template, named items are created in the
            output figure in addition to any items the figure
            already has in this array. You can modify these items
            in the output figure by making your own item with
            `templateitemname` matching this `name` alongside your
            modifications (including `visible: false` or `enabled:
            false` to hide it). Has no effect outside of a
            template.
        step
            The unit of measurement that the `count` value will set
            the range by.
        stepmode
            Sets the range update mode. If "backward", the range
            update shifts the start of range back "count" times
            "step" milliseconds. If "todate", the range update
            shifts the start of range back to the first timestamp
            from "count" times "step" milliseconds back. For
            example, with `step` set to "year" and `count` set to 1
            the range update shifts the start of the range back to
            January 01 of the current year. Month and year "todate"
            are currently available only for the built-in
            (Gregorian) calendar.
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
        visible
            Determines whether or not this button is visible.
        """

    def __init__(
        self,
        arg=None,
        count=None,
        label=None,
        name=None,
        step=None,
        stepmode=None,
        templateitemname=None,
        visible=None,
        **kwargs,
    ):
        """
        Construct a new Button object

        Sets the specifications for each buttons. By default, a range
        selector comes with no buttons.

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of :class:`plotly.graph_objs.layout.xaxis.r
            angeselector.Button`
        count
            Sets the number of steps to take to update the range.
            Use with `step` to specify the update interval.
        label
            Sets the text label to appear on the button.
        name
            When used in a template, named items are created in the
            output figure in addition to any items the figure
            already has in this array. You can modify these items
            in the output figure by making your own item with
            `templateitemname` matching this `name` alongside your
            modifications (including `visible: false` or `enabled:
            false` to hide it). Has no effect outside of a
            template.
        step
            The unit of measurement that the `count` value will set
            the range by.
        stepmode
            Sets the range update mode. If "backward", the range
            update shifts the start of range back "count" times
            "step" milliseconds. If "todate", the range update
            shifts the start of range back to the first timestamp
            from "count" times "step" milliseconds back. For
            example, with `step` set to "year" and `count` set to 1
            the range update shifts the start of the range back to
            January 01 of the current year. Month and year "todate"
            are currently available only for the built-in
            (Gregorian) calendar.
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
        visible
            Determines whether or not this button is visible.

        Returns
        -------
        Button
        """
        super(Button, self).__init__("buttons")

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
The first argument to the plotly.graph_objs.layout.xaxis.rangeselector.Button
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.xaxis.rangeselector.Button`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("count", None)
        _v = count if count is not None else _v
        if _v is not None:
            self["count"] = _v
        _v = arg.pop("label", None)
        _v = label if label is not None else _v
        if _v is not None:
            self["label"] = _v
        _v = arg.pop("name", None)
        _v = name if name is not None else _v
        if _v is not None:
            self["name"] = _v
        _v = arg.pop("step", None)
        _v = step if step is not None else _v
        if _v is not None:
            self["step"] = _v
        _v = arg.pop("stepmode", None)
        _v = stepmode if stepmode is not None else _v
        if _v is not None:
            self["stepmode"] = _v
        _v = arg.pop("templateitemname", None)
        _v = templateitemname if templateitemname is not None else _v
        if _v is not None:
            self["templateitemname"] = _v
        _v = arg.pop("visible", None)
        _v = visible if visible is not None else _v
        if _v is not None:
            self["visible"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
