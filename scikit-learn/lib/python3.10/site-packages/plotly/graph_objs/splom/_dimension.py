#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Dimension(_BaseTraceHierarchyType):
    _parent_path_str = "splom"
    _path_str = "splom.dimension"
    _valid_props = {
        "axis",
        "label",
        "name",
        "templateitemname",
        "values",
        "valuessrc",
        "visible",
    }

    @property
    def axis(self):
        """
        The 'axis' property is an instance of Axis
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.splom.dimension.Axis`
          - A dict of string/value properties that will be passed
            to the Axis constructor

        Returns
        -------
        plotly.graph_objs.splom.dimension.Axis
        """
        return self["axis"]

    @axis.setter
    def axis(self, val):
        self["axis"] = val

    @property
    def label(self):
        """
        Sets the label corresponding to this splom dimension.

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
    def values(self):
        """
        Sets the dimension values to be plotted.

        The 'values' property is an array that may be specified as a tuple,
        list, numpy array, or pandas Series

        Returns
        -------
        numpy.ndarray
        """
        return self["values"]

    @values.setter
    def values(self, val):
        self["values"] = val

    @property
    def valuessrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `values`.

        The 'valuessrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["valuessrc"]

    @valuessrc.setter
    def valuessrc(self, val):
        self["valuessrc"] = val

    @property
    def visible(self):
        """
        Determines whether or not this dimension is shown on the graph.
        Note that even visible false dimension contribute to the
        default grid generate by this splom trace.

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
        axis
            :class:`plotly.graph_objects.splom.dimension.Axis`
            instance or dict with compatible properties
        label
            Sets the label corresponding to this splom dimension.
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
        values
            Sets the dimension values to be plotted.
        valuessrc
            Sets the source reference on Chart Studio Cloud for
            `values`.
        visible
            Determines whether or not this dimension is shown on
            the graph. Note that even visible false dimension
            contribute to the default grid generate by this splom
            trace.
        """

    def __init__(
        self,
        arg=None,
        axis=None,
        label=None,
        name=None,
        templateitemname=None,
        values=None,
        valuessrc=None,
        visible=None,
        **kwargs,
    ):
        """
        Construct a new Dimension object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.splom.Dimension`
        axis
            :class:`plotly.graph_objects.splom.dimension.Axis`
            instance or dict with compatible properties
        label
            Sets the label corresponding to this splom dimension.
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
        values
            Sets the dimension values to be plotted.
        valuessrc
            Sets the source reference on Chart Studio Cloud for
            `values`.
        visible
            Determines whether or not this dimension is shown on
            the graph. Note that even visible false dimension
            contribute to the default grid generate by this splom
            trace.

        Returns
        -------
        Dimension
        """
        super().__init__("dimensions")
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
The first argument to the plotly.graph_objs.splom.Dimension
constructor must be a dict or
an instance of :class:`plotly.graph_objs.splom.Dimension`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("axis", arg, axis)
        self._set_property("label", arg, label)
        self._set_property("name", arg, name)
        self._set_property("templateitemname", arg, templateitemname)
        self._set_property("values", arg, values)
        self._set_property("valuessrc", arg, valuessrc)
        self._set_property("visible", arg, visible)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
