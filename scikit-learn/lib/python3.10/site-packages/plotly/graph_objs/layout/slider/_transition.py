#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Transition(_BaseLayoutHierarchyType):
    _parent_path_str = "layout.slider"
    _path_str = "layout.slider.transition"
    _valid_props = {"duration", "easing"}

    @property
    def duration(self):
        """
        Sets the duration of the slider transition

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

    @property
    def easing(self):
        """
        Sets the easing function of the slider transition

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

    @property
    def _prop_descriptions(self):
        return """\
        duration
            Sets the duration of the slider transition
        easing
            Sets the easing function of the slider transition
        """

    def __init__(self, arg=None, duration=None, easing=None, **kwargs):
        """
        Construct a new Transition object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.layout.slider.Transition`
        duration
            Sets the duration of the slider transition
        easing
            Sets the easing function of the slider transition

        Returns
        -------
        Transition
        """
        super().__init__("transition")
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
The first argument to the plotly.graph_objs.layout.slider.Transition
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.slider.Transition`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("duration", arg, duration)
        self._set_property("easing", arg, easing)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
