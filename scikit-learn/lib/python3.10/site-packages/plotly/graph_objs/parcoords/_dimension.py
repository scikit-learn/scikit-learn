#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Dimension(_BaseTraceHierarchyType):
    _parent_path_str = "parcoords"
    _path_str = "parcoords.dimension"
    _valid_props = {
        "constraintrange",
        "label",
        "multiselect",
        "name",
        "range",
        "templateitemname",
        "tickformat",
        "ticktext",
        "ticktextsrc",
        "tickvals",
        "tickvalssrc",
        "values",
        "valuessrc",
        "visible",
    }

    @property
    def constraintrange(self):
        """
            The domain range to which the filter on the dimension is
            constrained. Must be an array of `[fromValue, toValue]` with
            `fromValue <= toValue`, or if `multiselect` is not disabled,
            you may give an array of arrays, where each inner array is
            `[fromValue, toValue]`.

            The 'constraintrange' property is an info array that may be specified as:

            * a list or tuple of 2 elements where:
        (0) The 'constraintrange[0]' property accepts values of any type
        (1) The 'constraintrange[1]' property accepts values of any type

            * a 2D list where:
        (0) The 'constraintrange[i][0]' property accepts values of any type
        (1) The 'constraintrange[i][1]' property accepts values of any type

            Returns
            -------
            list
        """
        return self["constraintrange"]

    @constraintrange.setter
    def constraintrange(self, val):
        self["constraintrange"] = val

    @property
    def label(self):
        """
        The shown name of the dimension.

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
    def multiselect(self):
        """
        Do we allow multiple selection ranges or just a single range?

        The 'multiselect' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["multiselect"]

    @multiselect.setter
    def multiselect(self, val):
        self["multiselect"] = val

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
    def range(self):
        """
            The domain range that represents the full, shown axis extent.
            Defaults to the `values` extent. Must be an array of
            `[fromValue, toValue]` with finite numbers as elements.

            The 'range' property is an info array that may be specified as:

            * a list or tuple of 2 elements where:
        (0) The 'range[0]' property is a number and may be specified as:
              - An int or float
        (1) The 'range[1]' property is a number and may be specified as:
              - An int or float

            Returns
            -------
            list
        """
        return self["range"]

    @range.setter
    def range(self, val):
        self["range"] = val

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
    def tickformat(self):
        """
        Sets the tick label formatting rule using d3 formatting mini-
        languages which are very similar to those in Python. For
        numbers, see:
        https://github.com/d3/d3-format/tree/v1.4.5#d3-format. And for
        dates see: https://github.com/d3/d3-time-
        format/tree/v2.2.3#locale_format. We add two items to d3's date
        formatter: "%h" for half of the year as a decimal number as
        well as "%{n}f" for fractional seconds with n digits. For
        example, *2016-10-13 09:15:23.456* with tickformat
        "%H~%M~%S.%2f" would display "09~15~23.46"

        The 'tickformat' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["tickformat"]

    @tickformat.setter
    def tickformat(self, val):
        self["tickformat"] = val

    @property
    def ticktext(self):
        """
        Sets the text displayed at the ticks position via `tickvals`.

        The 'ticktext' property is an array that may be specified as a tuple,
        list, numpy array, or pandas Series

        Returns
        -------
        numpy.ndarray
        """
        return self["ticktext"]

    @ticktext.setter
    def ticktext(self, val):
        self["ticktext"] = val

    @property
    def ticktextsrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `ticktext`.

        The 'ticktextsrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["ticktextsrc"]

    @ticktextsrc.setter
    def ticktextsrc(self, val):
        self["ticktextsrc"] = val

    @property
    def tickvals(self):
        """
        Sets the values at which ticks on this axis appear.

        The 'tickvals' property is an array that may be specified as a tuple,
        list, numpy array, or pandas Series

        Returns
        -------
        numpy.ndarray
        """
        return self["tickvals"]

    @tickvals.setter
    def tickvals(self, val):
        self["tickvals"] = val

    @property
    def tickvalssrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `tickvals`.

        The 'tickvalssrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["tickvalssrc"]

    @tickvalssrc.setter
    def tickvalssrc(self, val):
        self["tickvalssrc"] = val

    @property
    def values(self):
        """
        Dimension values. `values[n]` represents the value of the `n`th
        point in the dataset, therefore the `values` vector for all
        dimensions must be the same (longer vectors will be truncated).
        Each value must be a finite number.

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
        Shows the dimension when set to `true` (the default). Hides the
        dimension for `false`.

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
        constraintrange
            The domain range to which the filter on the dimension
            is constrained. Must be an array of `[fromValue,
            toValue]` with `fromValue <= toValue`, or if
            `multiselect` is not disabled, you may give an array of
            arrays, where each inner array is `[fromValue,
            toValue]`.
        label
            The shown name of the dimension.
        multiselect
            Do we allow multiple selection ranges or just a single
            range?
        name
            When used in a template, named items are created in the
            output figure in addition to any items the figure
            already has in this array. You can modify these items
            in the output figure by making your own item with
            `templateitemname` matching this `name` alongside your
            modifications (including `visible: false` or `enabled:
            false` to hide it). Has no effect outside of a
            template.
        range
            The domain range that represents the full, shown axis
            extent. Defaults to the `values` extent. Must be an
            array of `[fromValue, toValue]` with finite numbers as
            elements.
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
        tickformat
            Sets the tick label formatting rule using d3 formatting
            mini-languages which are very similar to those in
            Python. For numbers, see:
            https://github.com/d3/d3-format/tree/v1.4.5#d3-format.
            And for dates see: https://github.com/d3/d3-time-
            format/tree/v2.2.3#locale_format. We add two items to
            d3's date formatter: "%h" for half of the year as a
            decimal number as well as "%{n}f" for fractional
            seconds with n digits. For example, *2016-10-13
            09:15:23.456* with tickformat "%H~%M~%S.%2f" would
            display "09~15~23.46"
        ticktext
            Sets the text displayed at the ticks position via
            `tickvals`.
        ticktextsrc
            Sets the source reference on Chart Studio Cloud for
            `ticktext`.
        tickvals
            Sets the values at which ticks on this axis appear.
        tickvalssrc
            Sets the source reference on Chart Studio Cloud for
            `tickvals`.
        values
            Dimension values. `values[n]` represents the value of
            the `n`th point in the dataset, therefore the `values`
            vector for all dimensions must be the same (longer
            vectors will be truncated). Each value must be a finite
            number.
        valuessrc
            Sets the source reference on Chart Studio Cloud for
            `values`.
        visible
            Shows the dimension when set to `true` (the default).
            Hides the dimension for `false`.
        """

    def __init__(
        self,
        arg=None,
        constraintrange=None,
        label=None,
        multiselect=None,
        name=None,
        range=None,
        templateitemname=None,
        tickformat=None,
        ticktext=None,
        ticktextsrc=None,
        tickvals=None,
        tickvalssrc=None,
        values=None,
        valuessrc=None,
        visible=None,
        **kwargs,
    ):
        """
        Construct a new Dimension object

        The dimensions (variables) of the parallel coordinates chart.
        2..60 dimensions are supported.

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.parcoords.Dimension`
        constraintrange
            The domain range to which the filter on the dimension
            is constrained. Must be an array of `[fromValue,
            toValue]` with `fromValue <= toValue`, or if
            `multiselect` is not disabled, you may give an array of
            arrays, where each inner array is `[fromValue,
            toValue]`.
        label
            The shown name of the dimension.
        multiselect
            Do we allow multiple selection ranges or just a single
            range?
        name
            When used in a template, named items are created in the
            output figure in addition to any items the figure
            already has in this array. You can modify these items
            in the output figure by making your own item with
            `templateitemname` matching this `name` alongside your
            modifications (including `visible: false` or `enabled:
            false` to hide it). Has no effect outside of a
            template.
        range
            The domain range that represents the full, shown axis
            extent. Defaults to the `values` extent. Must be an
            array of `[fromValue, toValue]` with finite numbers as
            elements.
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
        tickformat
            Sets the tick label formatting rule using d3 formatting
            mini-languages which are very similar to those in
            Python. For numbers, see:
            https://github.com/d3/d3-format/tree/v1.4.5#d3-format.
            And for dates see: https://github.com/d3/d3-time-
            format/tree/v2.2.3#locale_format. We add two items to
            d3's date formatter: "%h" for half of the year as a
            decimal number as well as "%{n}f" for fractional
            seconds with n digits. For example, *2016-10-13
            09:15:23.456* with tickformat "%H~%M~%S.%2f" would
            display "09~15~23.46"
        ticktext
            Sets the text displayed at the ticks position via
            `tickvals`.
        ticktextsrc
            Sets the source reference on Chart Studio Cloud for
            `ticktext`.
        tickvals
            Sets the values at which ticks on this axis appear.
        tickvalssrc
            Sets the source reference on Chart Studio Cloud for
            `tickvals`.
        values
            Dimension values. `values[n]` represents the value of
            the `n`th point in the dataset, therefore the `values`
            vector for all dimensions must be the same (longer
            vectors will be truncated). Each value must be a finite
            number.
        valuessrc
            Sets the source reference on Chart Studio Cloud for
            `values`.
        visible
            Shows the dimension when set to `true` (the default).
            Hides the dimension for `false`.

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
The first argument to the plotly.graph_objs.parcoords.Dimension
constructor must be a dict or
an instance of :class:`plotly.graph_objs.parcoords.Dimension`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("constraintrange", arg, constraintrange)
        self._set_property("label", arg, label)
        self._set_property("multiselect", arg, multiselect)
        self._set_property("name", arg, name)
        self._set_property("range", arg, range)
        self._set_property("templateitemname", arg, templateitemname)
        self._set_property("tickformat", arg, tickformat)
        self._set_property("ticktext", arg, ticktext)
        self._set_property("ticktextsrc", arg, ticktextsrc)
        self._set_property("tickvals", arg, tickvals)
        self._set_property("tickvalssrc", arg, tickvalssrc)
        self._set_property("values", arg, values)
        self._set_property("valuessrc", arg, valuessrc)
        self._set_property("visible", arg, visible)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
