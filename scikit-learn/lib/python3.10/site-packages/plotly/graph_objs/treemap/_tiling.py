#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Tiling(_BaseTraceHierarchyType):
    _parent_path_str = "treemap"
    _path_str = "treemap.tiling"
    _valid_props = {"flip", "packing", "pad", "squarifyratio"}

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

    @property
    def packing(self):
        """
        Determines d3 treemap solver. For more info please refer to
        https://github.com/d3/d3-hierarchy#treemap-tiling

        The 'packing' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['squarify', 'binary', 'dice', 'slice', 'slice-dice',
                'dice-slice']

        Returns
        -------
        Any
        """
        return self["packing"]

    @packing.setter
    def packing(self, val):
        self["packing"] = val

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

    @property
    def squarifyratio(self):
        """
        When using "squarify" `packing` algorithm, according to https:/
        /github.com/d3/d3-
        hierarchy/blob/v3.1.1/README.md#squarify_ratio this option
        specifies the desired aspect ratio of the generated rectangles.
        The ratio must be specified as a number greater than or equal
        to one. Note that the orientation of the generated rectangles
        (tall or wide) is not implied by the ratio; for example, a
        ratio of two will attempt to produce a mixture of rectangles
        whose width:height ratio is either 2:1 or 1:2. When using
        "squarify", unlike d3 which uses the Golden Ratio i.e.
        1.618034, Plotly applies 1 to increase squares in treemap
        layouts.

        The 'squarifyratio' property is a number and may be specified as:
          - An int or float in the interval [1, inf]

        Returns
        -------
        int|float
        """
        return self["squarifyratio"]

    @squarifyratio.setter
    def squarifyratio(self, val):
        self["squarifyratio"] = val

    @property
    def _prop_descriptions(self):
        return """\
        flip
            Determines if the positions obtained from solver are
            flipped on each axis.
        packing
            Determines d3 treemap solver. For more info please
            refer to https://github.com/d3/d3-hierarchy#treemap-
            tiling
        pad
            Sets the inner padding (in px).
        squarifyratio
            When using "squarify" `packing` algorithm, according to
            https://github.com/d3/d3-
            hierarchy/blob/v3.1.1/README.md#squarify_ratio this
            option specifies the desired aspect ratio of the
            generated rectangles. The ratio must be specified as a
            number greater than or equal to one. Note that the
            orientation of the generated rectangles (tall or wide)
            is not implied by the ratio; for example, a ratio of
            two will attempt to produce a mixture of rectangles
            whose width:height ratio is either 2:1 or 1:2. When
            using "squarify", unlike d3 which uses the Golden Ratio
            i.e. 1.618034, Plotly applies 1 to increase squares in
            treemap layouts.
        """

    def __init__(
        self, arg=None, flip=None, packing=None, pad=None, squarifyratio=None, **kwargs
    ):
        """
        Construct a new Tiling object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.treemap.Tiling`
        flip
            Determines if the positions obtained from solver are
            flipped on each axis.
        packing
            Determines d3 treemap solver. For more info please
            refer to https://github.com/d3/d3-hierarchy#treemap-
            tiling
        pad
            Sets the inner padding (in px).
        squarifyratio
            When using "squarify" `packing` algorithm, according to
            https://github.com/d3/d3-
            hierarchy/blob/v3.1.1/README.md#squarify_ratio this
            option specifies the desired aspect ratio of the
            generated rectangles. The ratio must be specified as a
            number greater than or equal to one. Note that the
            orientation of the generated rectangles (tall or wide)
            is not implied by the ratio; for example, a ratio of
            two will attempt to produce a mixture of rectangles
            whose width:height ratio is either 2:1 or 1:2. When
            using "squarify", unlike d3 which uses the Golden Ratio
            i.e. 1.618034, Plotly applies 1 to increase squares in
            treemap layouts.

        Returns
        -------
        Tiling
        """
        super().__init__("tiling")
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
The first argument to the plotly.graph_objs.treemap.Tiling
constructor must be a dict or
an instance of :class:`plotly.graph_objs.treemap.Tiling`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("flip", arg, flip)
        self._set_property("packing", arg, packing)
        self._set_property("pad", arg, pad)
        self._set_property("squarifyratio", arg, squarifyratio)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
