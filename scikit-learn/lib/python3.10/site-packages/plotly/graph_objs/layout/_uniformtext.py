#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Uniformtext(_BaseLayoutHierarchyType):
    _parent_path_str = "layout"
    _path_str = "layout.uniformtext"
    _valid_props = {"minsize", "mode"}

    @property
    def minsize(self):
        """
        Sets the minimum text size between traces of the same type.

        The 'minsize' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["minsize"]

    @minsize.setter
    def minsize(self, val):
        self["minsize"] = val

    @property
    def mode(self):
        """
        Determines how the font size for various text elements are
        uniformed between each trace type. If the computed text sizes
        were smaller than the minimum size defined by
        `uniformtext.minsize` using "hide" option hides the text; and
        using "show" option shows the text without further downscaling.
        Please note that if the size defined by `minsize` is greater
        than the font size defined by trace, then the `minsize` is
        used.

        The 'mode' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                [False, 'hide', 'show']

        Returns
        -------
        Any
        """
        return self["mode"]

    @mode.setter
    def mode(self, val):
        self["mode"] = val

    @property
    def _prop_descriptions(self):
        return """\
        minsize
            Sets the minimum text size between traces of the same
            type.
        mode
            Determines how the font size for various text elements
            are uniformed between each trace type. If the computed
            text sizes were smaller than the minimum size defined
            by `uniformtext.minsize` using "hide" option hides the
            text; and using "show" option shows the text without
            further downscaling. Please note that if the size
            defined by `minsize` is greater than the font size
            defined by trace, then the `minsize` is used.
        """

    def __init__(self, arg=None, minsize=None, mode=None, **kwargs):
        """
        Construct a new Uniformtext object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.layout.Uniformtext`
        minsize
            Sets the minimum text size between traces of the same
            type.
        mode
            Determines how the font size for various text elements
            are uniformed between each trace type. If the computed
            text sizes were smaller than the minimum size defined
            by `uniformtext.minsize` using "hide" option hides the
            text; and using "show" option shows the text without
            further downscaling. Please note that if the size
            defined by `minsize` is greater than the font size
            defined by trace, then the `minsize` is used.

        Returns
        -------
        Uniformtext
        """
        super().__init__("uniformtext")
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
The first argument to the plotly.graph_objs.layout.Uniformtext
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.Uniformtext`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("minsize", arg, minsize)
        self._set_property("mode", arg, mode)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
