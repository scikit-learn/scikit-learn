from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Button(_BaseLayoutHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "layout.updatemenu"
    _path_str = "layout.updatemenu.button"
    _valid_props = {
        "args",
        "args2",
        "execute",
        "label",
        "method",
        "name",
        "templateitemname",
        "visible",
    }

    # args
    # ----
    @property
    def args(self):
        """
            Sets the arguments values to be passed to the Plotly method set
            in `method` on click.

            The 'args' property is an info array that may be specified as:

            * a list or tuple of up to 3 elements where:
        (0) The 'args[0]' property accepts values of any type
        (1) The 'args[1]' property accepts values of any type
        (2) The 'args[2]' property accepts values of any type

            Returns
            -------
            list
        """
        return self["args"]

    @args.setter
    def args(self, val):
        self["args"] = val

    # args2
    # -----
    @property
    def args2(self):
        """
            Sets a 2nd set of `args`, these arguments values are passed to
            the Plotly method set in `method` when clicking this button
            while in the active state. Use this to create toggle buttons.

            The 'args2' property is an info array that may be specified as:

            * a list or tuple of up to 3 elements where:
        (0) The 'args2[0]' property accepts values of any type
        (1) The 'args2[1]' property accepts values of any type
        (2) The 'args2[2]' property accepts values of any type

            Returns
            -------
            list
        """
        return self["args2"]

    @args2.setter
    def args2(self, val):
        self["args2"] = val

    # execute
    # -------
    @property
    def execute(self):
        """
        When true, the API method is executed. When false, all other
        behaviors are the same and command execution is skipped. This
        may be useful when hooking into, for example, the
        `plotly_buttonclicked` method and executing the API command
        manually without losing the benefit of the updatemenu
        automatically binding to the state of the plot through the
        specification of `method` and `args`.

        The 'execute' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["execute"]

    @execute.setter
    def execute(self, val):
        self["execute"] = val

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

    # method
    # ------
    @property
    def method(self):
        """
        Sets the Plotly method to be called on click. If the `skip`
        method is used, the API updatemenu will function as normal but
        will perform no API calls and will not bind automatically to
        state updates. This may be used to create a component interface
        and attach to updatemenu events manually via JavaScript.

        The 'method' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['restyle', 'relayout', 'animate', 'update', 'skip']

        Returns
        -------
        Any
        """
        return self["method"]

    @method.setter
    def method(self, val):
        self["method"] = val

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
        args
            Sets the arguments values to be passed to the Plotly
            method set in `method` on click.
        args2
            Sets a 2nd set of `args`, these arguments values are
            passed to the Plotly method set in `method` when
            clicking this button while in the active state. Use
            this to create toggle buttons.
        execute
            When true, the API method is executed. When false, all
            other behaviors are the same and command execution is
            skipped. This may be useful when hooking into, for
            example, the `plotly_buttonclicked` method and
            executing the API command manually without losing the
            benefit of the updatemenu automatically binding to the
            state of the plot through the specification of `method`
            and `args`.
        label
            Sets the text label to appear on the button.
        method
            Sets the Plotly method to be called on click. If the
            `skip` method is used, the API updatemenu will function
            as normal but will perform no API calls and will not
            bind automatically to state updates. This may be used
            to create a component interface and attach to
            updatemenu events manually via JavaScript.
        name
            When used in a template, named items are created in the
            output figure in addition to any items the figure
            already has in this array. You can modify these items
            in the output figure by making your own item with
            `templateitemname` matching this `name` alongside your
            modifications (including `visible: false` or `enabled:
            false` to hide it). Has no effect outside of a
            template.
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
        args=None,
        args2=None,
        execute=None,
        label=None,
        method=None,
        name=None,
        templateitemname=None,
        visible=None,
        **kwargs,
    ):
        """
        Construct a new Button object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.layout.updatemenu.Button`
        args
            Sets the arguments values to be passed to the Plotly
            method set in `method` on click.
        args2
            Sets a 2nd set of `args`, these arguments values are
            passed to the Plotly method set in `method` when
            clicking this button while in the active state. Use
            this to create toggle buttons.
        execute
            When true, the API method is executed. When false, all
            other behaviors are the same and command execution is
            skipped. This may be useful when hooking into, for
            example, the `plotly_buttonclicked` method and
            executing the API command manually without losing the
            benefit of the updatemenu automatically binding to the
            state of the plot through the specification of `method`
            and `args`.
        label
            Sets the text label to appear on the button.
        method
            Sets the Plotly method to be called on click. If the
            `skip` method is used, the API updatemenu will function
            as normal but will perform no API calls and will not
            bind automatically to state updates. This may be used
            to create a component interface and attach to
            updatemenu events manually via JavaScript.
        name
            When used in a template, named items are created in the
            output figure in addition to any items the figure
            already has in this array. You can modify these items
            in the output figure by making your own item with
            `templateitemname` matching this `name` alongside your
            modifications (including `visible: false` or `enabled:
            false` to hide it). Has no effect outside of a
            template.
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
The first argument to the plotly.graph_objs.layout.updatemenu.Button
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.updatemenu.Button`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("args", None)
        _v = args if args is not None else _v
        if _v is not None:
            self["args"] = _v
        _v = arg.pop("args2", None)
        _v = args2 if args2 is not None else _v
        if _v is not None:
            self["args2"] = _v
        _v = arg.pop("execute", None)
        _v = execute if execute is not None else _v
        if _v is not None:
            self["execute"] = _v
        _v = arg.pop("label", None)
        _v = label if label is not None else _v
        if _v is not None:
            self["label"] = _v
        _v = arg.pop("method", None)
        _v = method if method is not None else _v
        if _v is not None:
            self["method"] = _v
        _v = arg.pop("name", None)
        _v = name if name is not None else _v
        if _v is not None:
            self["name"] = _v
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
