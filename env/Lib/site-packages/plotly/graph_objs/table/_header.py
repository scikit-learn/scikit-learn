from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Header(_BaseTraceHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "table"
    _path_str = "table.header"
    _valid_props = {
        "align",
        "alignsrc",
        "fill",
        "font",
        "format",
        "formatsrc",
        "height",
        "line",
        "prefix",
        "prefixsrc",
        "suffix",
        "suffixsrc",
        "values",
        "valuessrc",
    }

    # align
    # -----
    @property
    def align(self):
        """
        Sets the horizontal alignment of the `text` within the box. Has
        an effect only if `text` spans two or more lines (i.e. `text`
        contains one or more <br> HTML tags) or if an explicit width is
        set to override the text width.

        The 'align' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['left', 'center', 'right']
          - A tuple, list, or one-dimensional numpy array of the above

        Returns
        -------
        Any|numpy.ndarray
        """
        return self["align"]

    @align.setter
    def align(self, val):
        self["align"] = val

    # alignsrc
    # --------
    @property
    def alignsrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `align`.

        The 'alignsrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["alignsrc"]

    @alignsrc.setter
    def alignsrc(self, val):
        self["alignsrc"] = val

    # fill
    # ----
    @property
    def fill(self):
        """
        The 'fill' property is an instance of Fill
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.table.header.Fill`
          - A dict of string/value properties that will be passed
            to the Fill constructor

            Supported dict properties:

                color
                    Sets the cell fill color. It accepts either a
                    specific color or an array of colors or a 2D
                    array of colors.
                colorsrc
                    Sets the source reference on Chart Studio Cloud
                    for `color`.

        Returns
        -------
        plotly.graph_objs.table.header.Fill
        """
        return self["fill"]

    @fill.setter
    def fill(self, val):
        self["fill"] = val

    # font
    # ----
    @property
    def font(self):
        """
        The 'font' property is an instance of Font
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.table.header.Font`
          - A dict of string/value properties that will be passed
            to the Font constructor

            Supported dict properties:

                color

                colorsrc
                    Sets the source reference on Chart Studio Cloud
                    for `color`.
                family
                    HTML font family - the typeface that will be
                    applied by the web browser. The web browser
                    will only be able to apply a font if it is
                    available on the system which it operates.
                    Provide multiple font families, separated by
                    commas, to indicate the preference in which to
                    apply fonts if they aren't available on the
                    system. The Chart Studio Cloud (at
                    https://chart-studio.plotly.com or on-premise)
                    generates images on a server, where only a
                    select number of fonts are installed and
                    supported. These include "Arial", "Balto",
                    "Courier New", "Droid Sans", "Droid Serif",
                    "Droid Sans Mono", "Gravitas One", "Old
                    Standard TT", "Open Sans", "Overpass", "PT Sans
                    Narrow", "Raleway", "Times New Roman".
                familysrc
                    Sets the source reference on Chart Studio Cloud
                    for `family`.
                lineposition
                    Sets the kind of decoration line(s) with text,
                    such as an "under", "over" or "through" as well
                    as combinations e.g. "under+over", etc.
                linepositionsrc
                    Sets the source reference on Chart Studio Cloud
                    for `lineposition`.
                shadow
                    Sets the shape and color of the shadow behind
                    text. "auto" places minimal shadow and applies
                    contrast text font color. See
                    https://developer.mozilla.org/en-
                    US/docs/Web/CSS/text-shadow for additional
                    options.
                shadowsrc
                    Sets the source reference on Chart Studio Cloud
                    for `shadow`.
                size

                sizesrc
                    Sets the source reference on Chart Studio Cloud
                    for `size`.
                style
                    Sets whether a font should be styled with a
                    normal or italic face from its family.
                stylesrc
                    Sets the source reference on Chart Studio Cloud
                    for `style`.
                textcase
                    Sets capitalization of text. It can be used to
                    make text appear in all-uppercase or all-
                    lowercase, or with each word capitalized.
                textcasesrc
                    Sets the source reference on Chart Studio Cloud
                    for `textcase`.
                variant
                    Sets the variant of the font.
                variantsrc
                    Sets the source reference on Chart Studio Cloud
                    for `variant`.
                weight
                    Sets the weight (or boldness) of the font.
                weightsrc
                    Sets the source reference on Chart Studio Cloud
                    for `weight`.

        Returns
        -------
        plotly.graph_objs.table.header.Font
        """
        return self["font"]

    @font.setter
    def font(self, val):
        self["font"] = val

    # format
    # ------
    @property
    def format(self):
        """
        Sets the cell value formatting rule using d3 formatting mini-
        languages which are very similar to those in Python. For
        numbers, see:
        https://github.com/d3/d3-format/tree/v1.4.5#d3-format.

        The 'format' property is an array that may be specified as a tuple,
        list, numpy array, or pandas Series

        Returns
        -------
        numpy.ndarray
        """
        return self["format"]

    @format.setter
    def format(self, val):
        self["format"] = val

    # formatsrc
    # ---------
    @property
    def formatsrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `format`.

        The 'formatsrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["formatsrc"]

    @formatsrc.setter
    def formatsrc(self, val):
        self["formatsrc"] = val

    # height
    # ------
    @property
    def height(self):
        """
        The height of cells.

        The 'height' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["height"]

    @height.setter
    def height(self, val):
        self["height"] = val

    # line
    # ----
    @property
    def line(self):
        """
        The 'line' property is an instance of Line
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.table.header.Line`
          - A dict of string/value properties that will be passed
            to the Line constructor

            Supported dict properties:

                color

                colorsrc
                    Sets the source reference on Chart Studio Cloud
                    for `color`.
                width

                widthsrc
                    Sets the source reference on Chart Studio Cloud
                    for `width`.

        Returns
        -------
        plotly.graph_objs.table.header.Line
        """
        return self["line"]

    @line.setter
    def line(self, val):
        self["line"] = val

    # prefix
    # ------
    @property
    def prefix(self):
        """
        Prefix for cell values.

        The 'prefix' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string
          - A tuple, list, or one-dimensional numpy array of the above

        Returns
        -------
        str|numpy.ndarray
        """
        return self["prefix"]

    @prefix.setter
    def prefix(self, val):
        self["prefix"] = val

    # prefixsrc
    # ---------
    @property
    def prefixsrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `prefix`.

        The 'prefixsrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["prefixsrc"]

    @prefixsrc.setter
    def prefixsrc(self, val):
        self["prefixsrc"] = val

    # suffix
    # ------
    @property
    def suffix(self):
        """
        Suffix for cell values.

        The 'suffix' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string
          - A tuple, list, or one-dimensional numpy array of the above

        Returns
        -------
        str|numpy.ndarray
        """
        return self["suffix"]

    @suffix.setter
    def suffix(self, val):
        self["suffix"] = val

    # suffixsrc
    # ---------
    @property
    def suffixsrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `suffix`.

        The 'suffixsrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["suffixsrc"]

    @suffixsrc.setter
    def suffixsrc(self, val):
        self["suffixsrc"] = val

    # values
    # ------
    @property
    def values(self):
        """
        Header cell values. `values[m][n]` represents the value of the
        `n`th point in column `m`, therefore the `values[m]` vector
        length for all columns must be the same (longer vectors will be
        truncated). Each value must be a finite number or a string.

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

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        align
            Sets the horizontal alignment of the `text` within the
            box. Has an effect only if `text` spans two or more
            lines (i.e. `text` contains one or more <br> HTML tags)
            or if an explicit width is set to override the text
            width.
        alignsrc
            Sets the source reference on Chart Studio Cloud for
            `align`.
        fill
            :class:`plotly.graph_objects.table.header.Fill`
            instance or dict with compatible properties
        font
            :class:`plotly.graph_objects.table.header.Font`
            instance or dict with compatible properties
        format
            Sets the cell value formatting rule using d3 formatting
            mini-languages which are very similar to those in
            Python. For numbers, see:
            https://github.com/d3/d3-format/tree/v1.4.5#d3-format.
        formatsrc
            Sets the source reference on Chart Studio Cloud for
            `format`.
        height
            The height of cells.
        line
            :class:`plotly.graph_objects.table.header.Line`
            instance or dict with compatible properties
        prefix
            Prefix for cell values.
        prefixsrc
            Sets the source reference on Chart Studio Cloud for
            `prefix`.
        suffix
            Suffix for cell values.
        suffixsrc
            Sets the source reference on Chart Studio Cloud for
            `suffix`.
        values
            Header cell values. `values[m][n]` represents the value
            of the `n`th point in column `m`, therefore the
            `values[m]` vector length for all columns must be the
            same (longer vectors will be truncated). Each value
            must be a finite number or a string.
        valuessrc
            Sets the source reference on Chart Studio Cloud for
            `values`.
        """

    def __init__(
        self,
        arg=None,
        align=None,
        alignsrc=None,
        fill=None,
        font=None,
        format=None,
        formatsrc=None,
        height=None,
        line=None,
        prefix=None,
        prefixsrc=None,
        suffix=None,
        suffixsrc=None,
        values=None,
        valuessrc=None,
        **kwargs,
    ):
        """
        Construct a new Header object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of :class:`plotly.graph_objs.table.Header`
        align
            Sets the horizontal alignment of the `text` within the
            box. Has an effect only if `text` spans two or more
            lines (i.e. `text` contains one or more <br> HTML tags)
            or if an explicit width is set to override the text
            width.
        alignsrc
            Sets the source reference on Chart Studio Cloud for
            `align`.
        fill
            :class:`plotly.graph_objects.table.header.Fill`
            instance or dict with compatible properties
        font
            :class:`plotly.graph_objects.table.header.Font`
            instance or dict with compatible properties
        format
            Sets the cell value formatting rule using d3 formatting
            mini-languages which are very similar to those in
            Python. For numbers, see:
            https://github.com/d3/d3-format/tree/v1.4.5#d3-format.
        formatsrc
            Sets the source reference on Chart Studio Cloud for
            `format`.
        height
            The height of cells.
        line
            :class:`plotly.graph_objects.table.header.Line`
            instance or dict with compatible properties
        prefix
            Prefix for cell values.
        prefixsrc
            Sets the source reference on Chart Studio Cloud for
            `prefix`.
        suffix
            Suffix for cell values.
        suffixsrc
            Sets the source reference on Chart Studio Cloud for
            `suffix`.
        values
            Header cell values. `values[m][n]` represents the value
            of the `n`th point in column `m`, therefore the
            `values[m]` vector length for all columns must be the
            same (longer vectors will be truncated). Each value
            must be a finite number or a string.
        valuessrc
            Sets the source reference on Chart Studio Cloud for
            `values`.

        Returns
        -------
        Header
        """
        super(Header, self).__init__("header")

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
The first argument to the plotly.graph_objs.table.Header
constructor must be a dict or
an instance of :class:`plotly.graph_objs.table.Header`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("align", None)
        _v = align if align is not None else _v
        if _v is not None:
            self["align"] = _v
        _v = arg.pop("alignsrc", None)
        _v = alignsrc if alignsrc is not None else _v
        if _v is not None:
            self["alignsrc"] = _v
        _v = arg.pop("fill", None)
        _v = fill if fill is not None else _v
        if _v is not None:
            self["fill"] = _v
        _v = arg.pop("font", None)
        _v = font if font is not None else _v
        if _v is not None:
            self["font"] = _v
        _v = arg.pop("format", None)
        _v = format if format is not None else _v
        if _v is not None:
            self["format"] = _v
        _v = arg.pop("formatsrc", None)
        _v = formatsrc if formatsrc is not None else _v
        if _v is not None:
            self["formatsrc"] = _v
        _v = arg.pop("height", None)
        _v = height if height is not None else _v
        if _v is not None:
            self["height"] = _v
        _v = arg.pop("line", None)
        _v = line if line is not None else _v
        if _v is not None:
            self["line"] = _v
        _v = arg.pop("prefix", None)
        _v = prefix if prefix is not None else _v
        if _v is not None:
            self["prefix"] = _v
        _v = arg.pop("prefixsrc", None)
        _v = prefixsrc if prefixsrc is not None else _v
        if _v is not None:
            self["prefixsrc"] = _v
        _v = arg.pop("suffix", None)
        _v = suffix if suffix is not None else _v
        if _v is not None:
            self["suffix"] = _v
        _v = arg.pop("suffixsrc", None)
        _v = suffixsrc if suffixsrc is not None else _v
        if _v is not None:
            self["suffixsrc"] = _v
        _v = arg.pop("values", None)
        _v = values if values is not None else _v
        if _v is not None:
            self["values"] = _v
        _v = arg.pop("valuessrc", None)
        _v = valuessrc if valuessrc is not None else _v
        if _v is not None:
            self["valuessrc"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
