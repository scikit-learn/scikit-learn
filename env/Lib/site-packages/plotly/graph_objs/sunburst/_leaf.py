from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Leaf(_BaseTraceHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "sunburst"
    _path_str = "sunburst.leaf"
    _valid_props = {"opacity"}

    # opacity
    # -------
    @property
    def opacity(self):
        """
        Sets the opacity of the leaves. With colorscale it is defaulted
        to 1; otherwise it is defaulted to 0.7

        The 'opacity' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["opacity"]

    @opacity.setter
    def opacity(self, val):
        self["opacity"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        opacity
            Sets the opacity of the leaves. With colorscale it is
            defaulted to 1; otherwise it is defaulted to 0.7
        """

    def __init__(self, arg=None, opacity=None, **kwargs):
        """
        Construct a new Leaf object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of :class:`plotly.graph_objs.sunburst.Leaf`
        opacity
            Sets the opacity of the leaves. With colorscale it is
            defaulted to 1; otherwise it is defaulted to 0.7

        Returns
        -------
        Leaf
        """
        super(Leaf, self).__init__("leaf")

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
The first argument to the plotly.graph_objs.sunburst.Leaf
constructor must be a dict or
an instance of :class:`plotly.graph_objs.sunburst.Leaf`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("opacity", None)
        _v = opacity if opacity is not None else _v
        if _v is not None:
            self["opacity"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
