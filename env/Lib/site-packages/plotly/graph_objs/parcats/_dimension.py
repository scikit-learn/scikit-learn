from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Dimension(_BaseTraceHierarchyType):

    # class properties
    # --------------------
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

    # categoryarray
    # -------------
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

    # categoryarraysrc
    # ----------------
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

    # categoryorder
    # -------------
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

    # displayindex
    # ------------
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

    # label
    # -----
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

    # ticktext
    # --------
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

    # ticktextsrc
    # -----------
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

    # values
    # ------
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

    # valuessrc
    # ---------
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

    # visible
    # -------
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

    # Self properties description
    # ---------------------------
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
        super(Dimension, self).__init__("dimensions")

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
The first argument to the plotly.graph_objs.parcats.Dimension
constructor must be a dict or
an instance of :class:`plotly.graph_objs.parcats.Dimension`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("categoryarray", None)
        _v = categoryarray if categoryarray is not None else _v
        if _v is not None:
            self["categoryarray"] = _v
        _v = arg.pop("categoryarraysrc", None)
        _v = categoryarraysrc if categoryarraysrc is not None else _v
        if _v is not None:
            self["categoryarraysrc"] = _v
        _v = arg.pop("categoryorder", None)
        _v = categoryorder if categoryorder is not None else _v
        if _v is not None:
            self["categoryorder"] = _v
        _v = arg.pop("displayindex", None)
        _v = displayindex if displayindex is not None else _v
        if _v is not None:
            self["displayindex"] = _v
        _v = arg.pop("label", None)
        _v = label if label is not None else _v
        if _v is not None:
            self["label"] = _v
        _v = arg.pop("ticktext", None)
        _v = ticktext if ticktext is not None else _v
        if _v is not None:
            self["ticktext"] = _v
        _v = arg.pop("ticktextsrc", None)
        _v = ticktextsrc if ticktextsrc is not None else _v
        if _v is not None:
            self["ticktextsrc"] = _v
        _v = arg.pop("values", None)
        _v = values if values is not None else _v
        if _v is not None:
            self["values"] = _v
        _v = arg.pop("valuessrc", None)
        _v = valuessrc if valuessrc is not None else _v
        if _v is not None:
            self["valuessrc"] = _v
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
