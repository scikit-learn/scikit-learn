from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Transition(_BaseLayoutHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "layout"
    _path_str = "layout.transition"
    _valid_props = {"duration", "easing", "ordering"}

    # duration
    # --------
    @property
    def duration(self):
        """
        The duration of the transition, in milliseconds. If equal to
        zero, updates are synchronous.

        The 'duration' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["duration"]

    @duration.setter
    def duration(self, val):
        self["duration"] = val

    # easing
    # ------
    @property
    def easing(self):
        """
        The easing function used for the transition

        The 'easing' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['linear', 'quad', 'cubic', 'sin', 'exp', 'circle',
                'elastic', 'back', 'bounce', 'linear-in', 'quad-in',
                'cubic-in', 'sin-in', 'exp-in', 'circle-in', 'elastic-in',
                'back-in', 'bounce-in', 'linear-out', 'quad-out',
                'cubic-out', 'sin-out', 'exp-out', 'circle-out',
                'elastic-out', 'back-out', 'bounce-out', 'linear-in-out',
                'quad-in-out', 'cubic-in-out', 'sin-in-out', 'exp-in-out',
                'circle-in-out', 'elastic-in-out', 'back-in-out',
                'bounce-in-out']

        Returns
        -------
        Any
        """
        return self["easing"]

    @easing.setter
    def easing(self, val):
        self["easing"] = val

    # ordering
    # --------
    @property
    def ordering(self):
        """
        Determines whether the figure's layout or traces smoothly
        transitions during updates that make both traces and layout
        change.

        The 'ordering' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['layout first', 'traces first']

        Returns
        -------
        Any
        """
        return self["ordering"]

    @ordering.setter
    def ordering(self, val):
        self["ordering"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        duration
            The duration of the transition, in milliseconds. If
            equal to zero, updates are synchronous.
        easing
            The easing function used for the transition
        ordering
            Determines whether the figure's layout or traces
            smoothly transitions during updates that make both
            traces and layout change.
        """

    def __init__(self, arg=None, duration=None, easing=None, ordering=None, **kwargs):
        """
        Construct a new Transition object

        Sets transition options used during Plotly.react updates.

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.layout.Transition`
        duration
            The duration of the transition, in milliseconds. If
            equal to zero, updates are synchronous.
        easing
            The easing function used for the transition
        ordering
            Determines whether the figure's layout or traces
            smoothly transitions during updates that make both
            traces and layout change.

        Returns
        -------
        Transition
        """
        super(Transition, self).__init__("transition")

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
The first argument to the plotly.graph_objs.layout.Transition
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.Transition`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("duration", None)
        _v = duration if duration is not None else _v
        if _v is not None:
            self["duration"] = _v
        _v = arg.pop("easing", None)
        _v = easing if easing is not None else _v
        if _v is not None:
            self["easing"] = _v
        _v = arg.pop("ordering", None)
        _v = ordering if ordering is not None else _v
        if _v is not None:
            self["ordering"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
