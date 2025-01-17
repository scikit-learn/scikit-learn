from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Unselected(_BaseTraceHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "barpolar"
    _path_str = "barpolar.unselected"
    _valid_props = {"marker", "textfont"}

    # marker
    # ------
    @property
    def marker(self):
        """
        The 'marker' property is an instance of Marker
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.barpolar.unselected.Marker`
          - A dict of string/value properties that will be passed
            to the Marker constructor

            Supported dict properties:

                color
                    Sets the marker color of unselected points,
                    applied only when a selection exists.
                opacity
                    Sets the marker opacity of unselected points,
                    applied only when a selection exists.

        Returns
        -------
        plotly.graph_objs.barpolar.unselected.Marker
        """
        return self["marker"]

    @marker.setter
    def marker(self, val):
        self["marker"] = val

    # textfont
    # --------
    @property
    def textfont(self):
        """
        The 'textfont' property is an instance of Textfont
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.barpolar.unselected.Textfont`
          - A dict of string/value properties that will be passed
            to the Textfont constructor

            Supported dict properties:

                color
                    Sets the text font color of unselected points,
                    applied only when a selection exists.

        Returns
        -------
        plotly.graph_objs.barpolar.unselected.Textfont
        """
        return self["textfont"]

    @textfont.setter
    def textfont(self, val):
        self["textfont"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        marker
            :class:`plotly.graph_objects.barpolar.unselected.Marker
            ` instance or dict with compatible properties
        textfont
            :class:`plotly.graph_objects.barpolar.unselected.Textfo
            nt` instance or dict with compatible properties
        """

    def __init__(self, arg=None, marker=None, textfont=None, **kwargs):
        """
        Construct a new Unselected object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.barpolar.Unselected`
        marker
            :class:`plotly.graph_objects.barpolar.unselected.Marker
            ` instance or dict with compatible properties
        textfont
            :class:`plotly.graph_objects.barpolar.unselected.Textfo
            nt` instance or dict with compatible properties

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
The first argument to the plotly.graph_objs.barpolar.Unselected
constructor must be a dict or
an instance of :class:`plotly.graph_objs.barpolar.Unselected`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("marker", None)
        _v = marker if marker is not None else _v
        if _v is not None:
            self["marker"] = _v
        _v = arg.pop("textfont", None)
        _v = textfont if textfont is not None else _v
        if _v is not None:
            self["textfont"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
