#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Dimension(_BaseTraceHierarchyType):
    _parent_path_str = "parcats"
    _path_str = "parcats.dimension"
    _valid_props = {
        "categoryarray",
        "categoryarraysrc",
        "categoryorder",
        "displayindex",
        "label",
        "ticktext",
        "ticktextsrc",
        "values",
        "valuessrc",
        "visible",
    }

    @property
    def categoryarray(self):
        """
        Sets the order in which categories in this dimension appear.
        Only has an effect if `categoryorder` is set to "array". Used
        with `categoryorder`.

        The 'categoryarray' property is an array that may be specified as a tuple,
        list, numpy array, or pandas Series

        Returns
        -------
        numpy.ndarray
        """
        return self["categoryarray"]

    @categoryarray.setter
    def categoryarray(self, val):
        self["categoryarray"] = val

    @property
    def categoryarraysrc(self):
        """
        Sets the source reference on Chart Studio Cloud for
        `categoryarray`.

        The 'categoryarraysrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["categoryarraysrc"]

    @categoryarraysrc.setter
    def categoryarraysrc(self, val):
        self["categoryarraysrc"] = val

    @property
    def categoryorder(self):
        """
        Specifies the ordering logic for the categories in the
        dimension. By default, plotly uses "trace", which specifies the
        order that is present in the data supplied. Set `categoryorder`
        to *category ascending* or *category descending* if order
        should be determined by the alphanumerical order of the
        category names. Set `categoryorder` to "array" to derive the
        ordering from the attribute `categoryarray`. If a category is
        not found in the `categoryarray` array, the sorting behavior
        for that attribute will be identical to the "trace" mode. The
        unspecified categories will follow the categories in
        `categoryarray`.

        The 'categoryorder' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['trace', 'category ascending', 'category descending',
                'array']

        Returns
        -------
        Any
        """
        return self["categoryorder"]

    @categoryorder.setter
    def categoryorder(self, val):
        self["categoryorder"] = val

    @property
    def displayindex(self):
        """
        The display index of dimension, from left to right, zero
        indexed, defaults to dimension index.

        The 'displayindex' property is a integer and may be specified as:
          - An int (or float that will be cast to an int)

        Returns
        -------
        int
        """
        return self["displayindex"]

    @displayindex.setter
    def displayindex(self, val):
        self["displayindex"] = val

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
    def ticktext(self):
        """
        Sets alternative tick labels for the categories in this
        dimension. Only has an effect if `categoryorder` is set to
        "array". Should be an array the same length as `categoryarray`
        Used with `categoryorder`.

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
    def values(self):
        """
        Dimension values. `values[n]` represents the category value of
        the `n`th point in the dataset, therefore the `values` vector
        for all dimensions must be the same (longer vectors will be
        truncated).

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
        categoryarray
            Sets the order in which categories in this dimension
            appear. Only has an effect if `categoryorder` is set to
            "array". Used with `categoryorder`.
        categoryarraysrc
            Sets the source reference on Chart Studio Cloud for
            `categoryarray`.
        categoryorder
            Specifies the ordering logic for the categories in the
            dimension. By default, plotly uses "trace", which
            specifies the order that is present in the data
            supplied. Set `categoryorder` to *category ascending*
            or *category descending* if order should be determined
            by the alphanumerical order of the category names. Set
            `categoryorder` to "array" to derive the ordering from
            the attribute `categoryarray`. If a category is not
            found in the `categoryarray` array, the sorting
            behavior for that attribute will be identical to the
            "trace" mode. The unspecified categories will follow
            the categories in `categoryarray`.
        displayindex
            The display index of dimension, from left to right,
            zero indexed, defaults to dimension index.
        label
            The shown name of the dimension.
        ticktext
            Sets alternative tick labels for the categories in this
            dimension. Only has an effect if `categoryorder` is set
            to "array". Should be an array the same length as
            `categoryarray` Used with `categoryorder`.
        ticktextsrc
            Sets the source reference on Chart Studio Cloud for
            `ticktext`.
        values
            Dimension values. `values[n]` represents the category
            value of the `n`th point in the dataset, therefore the
            `values` vector for all dimensions must be the same
            (longer vectors will be truncated).
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
        categoryarray=None,
        categoryarraysrc=None,
        categoryorder=None,
        displayindex=None,
        label=None,
        ticktext=None,
        ticktextsrc=None,
        values=None,
        valuessrc=None,
        visible=None,
        **kwargs,
    ):
        """
        Construct a new Dimension object

        The dimensions (variables) of the parallel categories diagram.

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.parcats.Dimension`
        categoryarray
            Sets the order in which categories in this dimension
            appear. Only has an effect if `categoryorder` is set to
            "array". Used with `categoryorder`.
        categoryarraysrc
            Sets the source reference on Chart Studio Cloud for
            `categoryarray`.
        categoryorder
            Specifies the ordering logic for the categories in the
            dimension. By default, plotly uses "trace", which
            specifies the order that is present in the data
            supplied. Set `categoryorder` to *category ascending*
            or *category descending* if order should be determined
            by the alphanumerical order of the category names. Set
            `categoryorder` to "array" to derive the ordering from
            the attribute `categoryarray`. If a category is not
            found in the `categoryarray` array, the sorting
            behavior for that attribute will be identical to the
            "trace" mode. The unspecified categories will follow
            the categories in `categoryarray`.
        displayindex
            The display index of dimension, from left to right,
            zero indexed, defaults to dimension index.
        label
            The shown name of the dimension.
        ticktext
            Sets alternative tick labels for the categories in this
            dimension. Only has an effect if `categoryorder` is set
            to "array". Should be an array the same length as
            `categoryarray` Used with `categoryorder`.
        ticktextsrc
            Sets the source reference on Chart Studio Cloud for
            `ticktext`.
        values
            Dimension values. `values[n]` represents the category
            value of the `n`th point in the dataset, therefore the
            `values` vector for all dimensions must be the same
            (longer vectors will be truncated).
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
The first argument to the plotly.graph_objs.parcats.Dimension
constructor must be a dict or
an instance of :class:`plotly.graph_objs.parcats.Dimension`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("categoryarray", arg, categoryarray)
        self._set_property("categoryarraysrc", arg, categoryarraysrc)
        self._set_property("categoryorder", arg, categoryorder)
        self._set_property("displayindex", arg, displayindex)
        self._set_property("label", arg, label)
        self._set_property("ticktext", arg, ticktext)
        self._set_property("ticktextsrc", arg, ticktextsrc)
        self._set_property("values", arg, values)
        self._set_property("valuessrc", arg, valuessrc)
        self._set_property("visible", arg, visible)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
