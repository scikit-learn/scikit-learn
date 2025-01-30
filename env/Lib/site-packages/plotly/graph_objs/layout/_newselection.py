from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Newselection(_BaseLayoutHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "layout"
    _path_str = "layout.newselection"
    _valid_props = {"line", "mode"}

    # line
    # ----
    @property
    def line(self):
        """
        The 'line' property is an instance of Line
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.newselection.Line`
          - A dict of string/value properties that will be passed
            to the Line constructor

            Supported dict properties:

                color
                    Sets the line color. By default uses either
                    dark grey or white to increase contrast with
                    background color.
                dash
                    Sets the dash style of lines. Set to a dash
                    type string ("solid", "dot", "dash",
                    "longdash", "dashdot", or "longdashdot") or a
                    dash length list in px (eg "5px,10px,2px,2px").
                width
                    Sets the line width (in px).

        Returns
        -------
        plotly.graph_objs.layout.newselection.Line
        """
        return self["line"]

    @line.setter
    def line(self, val):
        self["line"] = val

    # mode
    # ----
    @property
    def mode(self):
        """
        Describes how a new selection is created. If `immediate`, a new
        selection is created after first mouse up. If `gradual`, a new
        selection is not created after first mouse. By adding to and
        subtracting from the initial selection, this option allows
        declaring extra outlines of the selection.

        The 'mode' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['immediate', 'gradual']

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
        line
            :class:`plotly.graph_objects.layout.newselection.Line`
            instance or dict with compatible properties
        mode
            Describes how a new selection is created. If
            `immediate`, a new selection is created after first
            mouse up. If `gradual`, a new selection is not created
            after first mouse. By adding to and subtracting from
            the initial selection, this option allows declaring
            extra outlines of the selection.
        """

    def __init__(self, arg=None, line=None, mode=None, **kwargs):
        """
        Construct a new Newselection object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.layout.Newselection`
        line
            :class:`plotly.graph_objects.layout.newselection.Line`
            instance or dict with compatible properties
        mode
            Describes how a new selection is created. If
            `immediate`, a new selection is created after first
            mouse up. If `gradual`, a new selection is not created
            after first mouse. By adding to and subtracting from
            the initial selection, this option allows declaring
            extra outlines of the selection.

        Returns
        -------
        Newselection
        """
        super(Newselection, self).__init__("newselection")

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
The first argument to the plotly.graph_objs.layout.Newselection
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.Newselection`"""
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
