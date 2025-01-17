from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Pathbar(_BaseTraceHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "icicle"
    _path_str = "icicle.pathbar"
    _valid_props = {"edgeshape", "side", "textfont", "thickness", "visible"}

    # edgeshape
    # ---------
    @property
    def edgeshape(self):
        """
        Determines which shape is used for edges between `barpath`
        labels.

        The 'edgeshape' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['>', '<', '|', '/', '\\']

        Returns
        -------
        Any
        """
        return self["edgeshape"]

    @edgeshape.setter
    def edgeshape(self, val):
        self["edgeshape"] = val

    # side
    # ----
    @property
    def side(self):
        """
        Determines on which side of the the treemap the `pathbar`
        should be presented.

        The 'side' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['top', 'bottom']

        Returns
        -------
        Any
        """
        return self["side"]

    @side.setter
    def side(self, val):
        self["side"] = val

    # textfont
    # --------
    @property
    def textfont(self):
        """
        Sets the font used inside `pathbar`.

        The 'textfont' property is an instance of Textfont
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.icicle.pathbar.Textfont`
          - A dict of string/value properties that will be passed
            to the Textfont constructor

            Supported dict properties:

                color

                colorsrc
                    Sets the source reference on Chart Studio Cloud
                    for `color`.
                family
                    HTML font family - the typeface that will be
                    applied by the web browser. The web browser
                    will only be able to apply a font if it is
                    available on the system which it operates.
                    Provide multiple font families, separated by
                    commas, to indicate the preference in which to
                    apply fonts if they aren't available on the
                    system. The Chart Studio Cloud (at
                    https://chart-studio.plotly.com or on-premise)
                    generates images on a server, where only a
                    select number of fonts are installed and
                    supported. These include "Arial", "Balto",
                    "Courier New", "Droid Sans", "Droid Serif",
                    "Droid Sans Mono", "Gravitas One", "Old
                    Standard TT", "Open Sans", "Overpass", "PT Sans
                    Narrow", "Raleway", "Times New Roman".
                familysrc
                    Sets the source reference on Chart Studio Cloud
                    for `family`.
                lineposition
                    Sets the kind of decoration line(s) with text,
                    such as an "under", "over" or "through" as well
                    as combinations e.g. "under+over", etc.
                linepositionsrc
                    Sets the source reference on Chart Studio Cloud
                    for `lineposition`.
                shadow
                    Sets the shape and color of the shadow behind
                    text. "auto" places minimal shadow and applies
                    contrast text font color. See
                    https://developer.mozilla.org/en-
                    US/docs/Web/CSS/text-shadow for additional
                    options.
                shadowsrc
                    Sets the source reference on Chart Studio Cloud
                    for `shadow`.
                size

                sizesrc
                    Sets the source reference on Chart Studio Cloud
                    for `size`.
                style
                    Sets whether a font should be styled with a
                    normal or italic face from its family.
                stylesrc
                    Sets the source reference on Chart Studio Cloud
                    for `style`.
                textcase
                    Sets capitalization of text. It can be used to
                    make text appear in all-uppercase or all-
                    lowercase, or with each word capitalized.
                textcasesrc
                    Sets the source reference on Chart Studio Cloud
                    for `textcase`.
                variant
                    Sets the variant of the font.
                variantsrc
                    Sets the source reference on Chart Studio Cloud
                    for `variant`.
                weight
                    Sets the weight (or boldness) of the font.
                weightsrc
                    Sets the source reference on Chart Studio Cloud
                    for `weight`.

        Returns
        -------
        plotly.graph_objs.icicle.pathbar.Textfont
        """
        return self["textfont"]

    @textfont.setter
    def textfont(self, val):
        self["textfont"] = val

    # thickness
    # ---------
    @property
    def thickness(self):
        """
        Sets the thickness of `pathbar` (in px). If not specified the
        `pathbar.textfont.size` is used with 3 pixles extra padding on
        each side.

        The 'thickness' property is a number and may be specified as:
          - An int or float in the interval [12, inf]

        Returns
        -------
        int|float
        """
        return self["thickness"]

    @thickness.setter
    def thickness(self, val):
        self["thickness"] = val

    # visible
    # -------
    @property
    def visible(self):
        """
        Determines if the path bar is drawn i.e. outside the trace
        `domain` and with one pixel gap.

        The 'visible' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["visible"]

    @visible.setter
    def visible(self, val):
        self["visible"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        edgeshape
            Determines which shape is used for edges between
            `barpath` labels.
        side
            Determines on which side of the the treemap the
            `pathbar` should be presented.
        textfont
            Sets the font used inside `pathbar`.
        thickness
            Sets the thickness of `pathbar` (in px). If not
            specified the `pathbar.textfont.size` is used with 3
            pixles extra padding on each side.
        visible
            Determines if the path bar is drawn i.e. outside the
            trace `domain` and with one pixel gap.
        """

    def __init__(
        self,
        arg=None,
        edgeshape=None,
        side=None,
        textfont=None,
        thickness=None,
        visible=None,
        **kwargs,
    ):
        """
        Construct a new Pathbar object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.icicle.Pathbar`
        edgeshape
            Determines which shape is used for edges between
            `barpath` labels.
        side
            Determines on which side of the the treemap the
            `pathbar` should be presented.
        textfont
            Sets the font used inside `pathbar`.
        thickness
            Sets the thickness of `pathbar` (in px). If not
            specified the `pathbar.textfont.size` is used with 3
            pixles extra padding on each side.
        visible
            Determines if the path bar is drawn i.e. outside the
            trace `domain` and with one pixel gap.

        Returns
        -------
        Pathbar
        """
        super(Pathbar, self).__init__("pathbar")

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
The first argument to the plotly.graph_objs.icicle.Pathbar
constructor must be a dict or
an instance of :class:`plotly.graph_objs.icicle.Pathbar`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("edgeshape", None)
        _v = edgeshape if edgeshape is not None else _v
        if _v is not None:
            self["edgeshape"] = _v
        _v = arg.pop("side", None)
        _v = side if side is not None else _v
        if _v is not None:
            self["side"] = _v
        _v = arg.pop("textfont", None)
        _v = textfont if textfont is not None else _v
        if _v is not None:
            self["textfont"] = _v
        _v = arg.pop("thickness", None)
        _v = thickness if thickness is not None else _v
        if _v is not None:
            self["thickness"] = _v
        _v = arg.pop("visible", None)
        _v = visible if visible is not None else _v
        if _v is not None:
            self["visible"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
