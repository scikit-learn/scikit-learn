from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Pad(_BaseLayoutHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "layout.updatemenu"
    _path_str = "layout.updatemenu.pad"
    _valid_props = {"b", "l", "r", "t"}

    # b
    # -
    @property
    def b(self):
        """
        The amount of padding (in px) along the bottom of the
        component.

        The 'b' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["b"]

    @b.setter
    def b(self, val):
        self["b"] = val

    # l
    # -
    @property
    def l(self):
        """
        The amount of padding (in px) on the left side of the
        component.

        The 'l' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["l"]

    @l.setter
    def l(self, val):
        self["l"] = val

    # r
    # -
    @property
    def r(self):
        """
        The amount of padding (in px) on the right side of the
        component.

        The 'r' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["r"]

    @r.setter
    def r(self, val):
        self["r"] = val

    # t
    # -
    @property
    def t(self):
        """
        The amount of padding (in px) along the top of the component.

        The 't' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["t"]

    @t.setter
    def t(self, val):
        self["t"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        b
            The amount of padding (in px) along the bottom of the
            component.
        l
            The amount of padding (in px) on the left side of the
            component.
        r
            The amount of padding (in px) on the right side of the
            component.
        t
            The amount of padding (in px) along the top of the
            component.
        """

    def __init__(self, arg=None, b=None, l=None, r=None, t=None, **kwargs):
        """
        Construct a new Pad object

        Sets the padding around the buttons or dropdown menu.

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.layout.updatemenu.Pad`
        b
            The amount of padding (in px) along the bottom of the
            component.
        l
            The amount of padding (in px) on the left side of the
            component.
        r
            The amount of padding (in px) on the right side of the
            component.
        t
            The amount of padding (in px) along the top of the
            component.

        Returns
        -------
        Pad
        """
        super(Pad, self).__init__("pad")

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
The first argument to the plotly.graph_objs.layout.updatemenu.Pad
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.updatemenu.Pad`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("b", None)
        _v = b if b is not None else _v
        if _v is not None:
            self["b"] = _v
        _v = arg.pop("l", None)
        _v = l if l is not None else _v
        if _v is not None:
            self["l"] = _v
        _v = arg.pop("r", None)
        _v = r if r is not None else _v
        if _v is not None:
            self["r"] = _v
        _v = arg.pop("t", None)
        _v = t if t is not None else _v
        if _v is not None:
            self["t"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
