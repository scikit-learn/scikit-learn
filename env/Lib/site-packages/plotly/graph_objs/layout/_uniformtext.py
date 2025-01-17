from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Uniformtext(_BaseLayoutHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "layout"
    _path_str = "layout.uniformtext"
    _valid_props = {"minsize", "mode"}

    # minsize
    # -------
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

    # mode
    # ----
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

    # Self properties description
    # ---------------------------
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
        super(Uniformtext, self).__init__("uniformtext")

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
The first argument to the plotly.graph_objs.layout.Uniformtext
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.Uniformtext`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("minsize", None)
        _v = minsize if minsize is not None else _v
        if _v is not None:
            self["minsize"] = _v
        _v = arg.pop("mode", None)
        _v = mode if mode is not None else _v
        if _v is not None:
            self["mode"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
