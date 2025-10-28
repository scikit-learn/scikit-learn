#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Spaceframe(_BaseTraceHierarchyType):
    _parent_path_str = "isosurface"
    _path_str = "isosurface.spaceframe"
    _valid_props = {"fill", "show"}

    @property
    def fill(self):
        """
        Sets the fill ratio of the `spaceframe` elements. The default
        fill value is 0.15 meaning that only 15% of the area of every
        faces of tetras would be shaded. Applying a greater `fill`
        ratio would allow the creation of stronger elements or could be
        sued to have entirely closed areas (in case of using 1).

        The 'fill' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["fill"]

    @fill.setter
    def fill(self, val):
        self["fill"] = val

    @property
    def show(self):
        """
        Displays/hides tetrahedron shapes between minimum and maximum
        iso-values. Often useful when either caps or surfaces are
        disabled or filled with values less than 1.

        The 'show' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["show"]

    @show.setter
    def show(self, val):
        self["show"] = val

    @property
    def _prop_descriptions(self):
        return """\
        fill
            Sets the fill ratio of the `spaceframe` elements. The
            default fill value is 0.15 meaning that only 15% of the
            area of every faces of tetras would be shaded. Applying
            a greater `fill` ratio would allow the creation of
            stronger elements or could be sued to have entirely
            closed areas (in case of using 1).
        show
            Displays/hides tetrahedron shapes between minimum and
            maximum iso-values. Often useful when either caps or
            surfaces are disabled or filled with values less than
            1.
        """

    def __init__(self, arg=None, fill=None, show=None, **kwargs):
        """
        Construct a new Spaceframe object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.isosurface.Spaceframe`
        fill
            Sets the fill ratio of the `spaceframe` elements. The
            default fill value is 0.15 meaning that only 15% of the
            area of every faces of tetras would be shaded. Applying
            a greater `fill` ratio would allow the creation of
            stronger elements or could be sued to have entirely
            closed areas (in case of using 1).
        show
            Displays/hides tetrahedron shapes between minimum and
            maximum iso-values. Often useful when either caps or
            surfaces are disabled or filled with values less than
            1.

        Returns
        -------
        Spaceframe
        """
        super().__init__("spaceframe")
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
The first argument to the plotly.graph_objs.isosurface.Spaceframe
constructor must be a dict or
an instance of :class:`plotly.graph_objs.isosurface.Spaceframe`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("fill", arg, fill)
        self._set_property("show", arg, show)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
