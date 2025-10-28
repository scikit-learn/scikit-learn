#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class X(_BaseTraceHierarchyType):
    _parent_path_str = "volume.caps"
    _path_str = "volume.caps.x"
    _valid_props = {"fill", "show"}

    @property
    def fill(self):
        """
        Sets the fill ratio of the `caps`. The default fill value of
        the `caps` is 1 meaning that they are entirely shaded. On the
        other hand Applying a `fill` ratio less than one would allow
        the creation of openings parallel to the edges.

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
        Sets the fill ratio of the `slices`. The default fill value of
        the x `slices` is 1 meaning that they are entirely shaded. On
        the other hand Applying a `fill` ratio less than one would
        allow the creation of openings parallel to the edges.

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
            Sets the fill ratio of the `caps`. The default fill
            value of the `caps` is 1 meaning that they are entirely
            shaded. On the other hand Applying a `fill` ratio less
            than one would allow the creation of openings parallel
            to the edges.
        show
            Sets the fill ratio of the `slices`. The default fill
            value of the x `slices` is 1 meaning that they are
            entirely shaded. On the other hand Applying a `fill`
            ratio less than one would allow the creation of
            openings parallel to the edges.
        """

    def __init__(self, arg=None, fill=None, show=None, **kwargs):
        """
        Construct a new X object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of :class:`plotly.graph_objs.volume.caps.X`
        fill
            Sets the fill ratio of the `caps`. The default fill
            value of the `caps` is 1 meaning that they are entirely
            shaded. On the other hand Applying a `fill` ratio less
            than one would allow the creation of openings parallel
            to the edges.
        show
            Sets the fill ratio of the `slices`. The default fill
            value of the x `slices` is 1 meaning that they are
            entirely shaded. On the other hand Applying a `fill`
            ratio less than one would allow the creation of
            openings parallel to the edges.

        Returns
        -------
        X
        """
        super().__init__("x")
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
The first argument to the plotly.graph_objs.volume.caps.X
constructor must be a dict or
an instance of :class:`plotly.graph_objs.volume.caps.X`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("fill", arg, fill)
        self._set_property("show", arg, show)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
