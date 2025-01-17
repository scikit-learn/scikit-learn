from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Tiling(_BaseTraceHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "icicle"
    _path_str = "icicle.tiling"
    _valid_props = {"flip", "orientation", "pad"}

    # flip
    # ----
    @property
    def flip(self):
        """
        Determines if the positions obtained from solver are flipped on
        each axis.

        The 'flip' property is a flaglist and may be specified
        as a string containing:
          - Any combination of ['x', 'y'] joined with '+' characters
            (e.g. 'x+y')

        Returns
        -------
        Any
        """
        return self["flip"]

    @flip.setter
    def flip(self, val):
        self["flip"] = val

    # orientation
    # -----------
    @property
    def orientation(self):
        """
        When set in conjunction with `tiling.flip`, determines on which
        side the root nodes are drawn in the chart. If
        `tiling.orientation` is "v" and `tiling.flip` is "", the root
        nodes appear at the top. If `tiling.orientation` is "v" and
        `tiling.flip` is "y", the root nodes appear at the bottom. If
        `tiling.orientation` is "h" and `tiling.flip` is "", the root
        nodes appear at the left. If `tiling.orientation` is "h" and
        `tiling.flip` is "x", the root nodes appear at the right.

        The 'orientation' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['v', 'h']

        Returns
        -------
        Any
        """
        return self["orientation"]

    @orientation.setter
    def orientation(self, val):
        self["orientation"] = val

    # pad
    # ---
    @property
    def pad(self):
        """
        Sets the inner padding (in px).

        The 'pad' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["pad"]

    @pad.setter
    def pad(self, val):
        self["pad"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        flip
            Determines if the positions obtained from solver are
            flipped on each axis.
        orientation
            When set in conjunction with `tiling.flip`, determines
            on which side the root nodes are drawn in the chart. If
            `tiling.orientation` is "v" and `tiling.flip` is "",
            the root nodes appear at the top. If
            `tiling.orientation` is "v" and `tiling.flip` is "y",
            the root nodes appear at the bottom. If
            `tiling.orientation` is "h" and `tiling.flip` is "",
            the root nodes appear at the left. If
            `tiling.orientation` is "h" and `tiling.flip` is "x",
            the root nodes appear at the right.
        pad
            Sets the inner padding (in px).
        """

    def __init__(self, arg=None, flip=None, orientation=None, pad=None, **kwargs):
        """
        Construct a new Tiling object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of :class:`plotly.graph_objs.icicle.Tiling`
        flip
            Determines if the positions obtained from solver are
            flipped on each axis.
        orientation
            When set in conjunction with `tiling.flip`, determines
            on which side the root nodes are drawn in the chart. If
            `tiling.orientation` is "v" and `tiling.flip` is "",
            the root nodes appear at the top. If
            `tiling.orientation` is "v" and `tiling.flip` is "y",
            the root nodes appear at the bottom. If
            `tiling.orientation` is "h" and `tiling.flip` is "",
            the root nodes appear at the left. If
            `tiling.orientation` is "h" and `tiling.flip` is "x",
            the root nodes appear at the right.
        pad
            Sets the inner padding (in px).

        Returns
        -------
        Tiling
        """
        super(Tiling, self).__init__("tiling")

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
The first argument to the plotly.graph_objs.icicle.Tiling
constructor must be a dict or
an instance of :class:`plotly.graph_objs.icicle.Tiling`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("flip", None)
        _v = flip if flip is not None else _v
        if _v is not None:
            self["flip"] = _v
        _v = arg.pop("orientation", None)
        _v = orientation if orientation is not None else _v
        if _v is not None:
            self["orientation"] = _v
        _v = arg.pop("pad", None)
        _v = pad if pad is not None else _v
        if _v is not None:
            self["pad"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
