from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Unselected(_BaseTraceHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "parcoords"
    _path_str = "parcoords.unselected"
    _valid_props = {"line"}

    # line
    # ----
    @property
    def line(self):
        """
        The 'line' property is an instance of Line
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.parcoords.unselected.Line`
          - A dict of string/value properties that will be passed
            to the Line constructor

            Supported dict properties:

                color
                    Sets the base color of unselected lines. in
                    connection with `unselected.line.opacity`.
                opacity
                    Sets the opacity of unselected lines. The
                    default "auto" decreases the opacity smoothly
                    as the number of lines increases. Use 1 to
                    achieve exact `unselected.line.color`.

        Returns
        -------
        plotly.graph_objs.parcoords.unselected.Line
        """
        return self["line"]

    @line.setter
    def line(self, val):
        self["line"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        line
            :class:`plotly.graph_objects.parcoords.unselected.Line`
            instance or dict with compatible properties
        """

    def __init__(self, arg=None, line=None, **kwargs):
        """
        Construct a new Unselected object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.parcoords.Unselected`
        line
            :class:`plotly.graph_objects.parcoords.unselected.Line`
            instance or dict with compatible properties

        Returns
        -------
        Unselected
        """
        super(Unselected, self).__init__("unselected")

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
The first argument to the plotly.graph_objs.parcoords.Unselected
constructor must be a dict or
an instance of :class:`plotly.graph_objs.parcoords.Unselected`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("line", None)
        _v = line if line is not None else _v
        if _v is not None:
            self["line"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
