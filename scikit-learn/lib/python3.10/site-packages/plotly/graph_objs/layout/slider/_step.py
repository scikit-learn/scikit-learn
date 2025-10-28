#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Step(_BaseLayoutHierarchyType):
    _parent_path_str = "layout.slider"
    _path_str = "layout.slider.step"
    _valid_props = {
        "args",
        "execute",
        "label",
        "method",
        "name",
        "templateitemname",
        "value",
        "visible",
    }

    @property
    def args(self):
        """
            Sets the arguments values to be passed to the Plotly method set
            in `method` on slide.

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

    @property
    def execute(self):
        """
        When true, the API method is executed. When false, all other
        behaviors are the same and command execution is skipped. This
        may be useful when hooking into, for example, the
        `plotly_sliderchange` method and executing the API command
        manually without losing the benefit of the slider automatically
        binding to the state of the plot through the specification of
        `method` and `args`.

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

    @property
    def label(self):
        """
        Sets the text label to appear on the slider

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

    @property
    def method(self):
        """
        Sets the Plotly method to be called when the slider value is
        changed. If the `skip` method is used, the API slider will
        function as normal but will perform no API calls and will not
        bind automatically to state updates. This may be used to create
        a component interface and attach to slider events manually via
        JavaScript.

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

    @property
    def value(self):
        """
        Sets the value of the slider step, used to refer to the step
        programatically. Defaults to the slider label if not provided.

        The 'value' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["value"]

    @value.setter
    def value(self, val):
        self["value"] = val

    @property
    def visible(self):
        """
        Determines whether or not this step is included in the slider.

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

    @property
    def _prop_descriptions(self):
        return """\
        args
            Sets the arguments values to be passed to the Plotly
            method set in `method` on slide.
        execute
            When true, the API method is executed. When false, all
            other behaviors are the same and command execution is
            skipped. This may be useful when hooking into, for
            example, the `plotly_sliderchange` method and executing
            the API command manually without losing the benefit of
            the slider automatically binding to the state of the
            plot through the specification of `method` and `args`.
        label
            Sets the text label to appear on the slider
        method
            Sets the Plotly method to be called when the slider
            value is changed. If the `skip` method is used, the API
            slider will function as normal but will perform no API
            calls and will not bind automatically to state updates.
            This may be used to create a component interface and
            attach to slider events manually via JavaScript.
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
        value
            Sets the value of the slider step, used to refer to the
            step programatically. Defaults to the slider label if
            not provided.
        visible
            Determines whether or not this step is included in the
            slider.
        """

    def __init__(
        self,
        arg=None,
        args=None,
        execute=None,
        label=None,
        method=None,
        name=None,
        templateitemname=None,
        value=None,
        visible=None,
        **kwargs,
    ):
        """
        Construct a new Step object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.layout.slider.Step`
        args
            Sets the arguments values to be passed to the Plotly
            method set in `method` on slide.
        execute
            When true, the API method is executed. When false, all
            other behaviors are the same and command execution is
            skipped. This may be useful when hooking into, for
            example, the `plotly_sliderchange` method and executing
            the API command manually without losing the benefit of
            the slider automatically binding to the state of the
            plot through the specification of `method` and `args`.
        label
            Sets the text label to appear on the slider
        method
            Sets the Plotly method to be called when the slider
            value is changed. If the `skip` method is used, the API
            slider will function as normal but will perform no API
            calls and will not bind automatically to state updates.
            This may be used to create a component interface and
            attach to slider events manually via JavaScript.
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
        value
            Sets the value of the slider step, used to refer to the
            step programatically. Defaults to the slider label if
            not provided.
        visible
            Determines whether or not this step is included in the
            slider.

        Returns
        -------
        Step
        """
        super().__init__("steps")
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
The first argument to the plotly.graph_objs.layout.slider.Step
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.slider.Step`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("args", arg, args)
        self._set_property("execute", arg, execute)
        self._set_property("label", arg, label)
        self._set_property("method", arg, method)
        self._set_property("name", arg, name)
        self._set_property("templateitemname", arg, templateitemname)
        self._set_property("value", arg, value)
        self._set_property("visible", arg, visible)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
