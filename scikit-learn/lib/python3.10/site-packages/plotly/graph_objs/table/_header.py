#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Header(_BaseTraceHierarchyType):
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

    @property
    def fill(self):
        """
        The 'fill' property is an instance of Fill
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.table.header.Fill`
          - A dict of string/value properties that will be passed
            to the Fill constructor

        Returns
        -------
        plotly.graph_objs.table.header.Fill
        """
        return self["fill"]

    @fill.setter
    def fill(self, val):
        self["fill"] = val

    @property
    def font(self):
        """
        The 'font' property is an instance of Font
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.table.header.Font`
          - A dict of string/value properties that will be passed
            to the Font constructor

        Returns
        -------
        plotly.graph_objs.table.header.Font
        """
        return self["font"]

    @font.setter
    def font(self, val):
        self["font"] = val

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

    @property
    def line(self):
        """
        The 'line' property is an instance of Line
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.table.header.Line`
          - A dict of string/value properties that will be passed
            to the Line constructor

        Returns
        -------
        plotly.graph_objs.table.header.Line
        """
        return self["line"]

    @line.setter
    def line(self, val):
        self["line"] = val

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
        super().__init__("header")
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
The first argument to the plotly.graph_objs.table.Header
constructor must be a dict or
an instance of :class:`plotly.graph_objs.table.Header`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("align", arg, align)
        self._set_property("alignsrc", arg, alignsrc)
        self._set_property("fill", arg, fill)
        self._set_property("font", arg, font)
        self._set_property("format", arg, format)
        self._set_property("formatsrc", arg, formatsrc)
        self._set_property("height", arg, height)
        self._set_property("line", arg, line)
        self._set_property("prefix", arg, prefix)
        self._set_property("prefixsrc", arg, prefixsrc)
        self._set_property("suffix", arg, suffix)
        self._set_property("suffixsrc", arg, suffixsrc)
        self._set_property("values", arg, values)
        self._set_property("valuessrc", arg, valuessrc)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
