from plotly.basedatatypes import BaseTraceHierarchyType as _BaseTraceHierarchyType
import copy as _copy


class Pattern(_BaseTraceHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "sunburst.marker"
    _path_str = "sunburst.marker.pattern"
    _valid_props = {
        "bgcolor",
        "bgcolorsrc",
        "fgcolor",
        "fgcolorsrc",
        "fgopacity",
        "fillmode",
        "shape",
        "shapesrc",
        "size",
        "sizesrc",
        "solidity",
        "soliditysrc",
    }

    # bgcolor
    # -------
    @property
    def bgcolor(self):
        """
        When there is no colorscale sets the color of background
        pattern fill. Defaults to a `marker.color` background when
        `fillmode` is "overlay". Otherwise, defaults to a transparent
        background.

        The 'bgcolor' property is a color and may be specified as:
          - A hex string (e.g. '#ff0000')
          - An rgb/rgba string (e.g. 'rgb(255,0,0)')
          - An hsl/hsla string (e.g. 'hsl(0,100%,50%)')
          - An hsv/hsva string (e.g. 'hsv(0,100%,100%)')
          - A named CSS color:
                aliceblue, antiquewhite, aqua, aquamarine, azure,
                beige, bisque, black, blanchedalmond, blue,
                blueviolet, brown, burlywood, cadetblue,
                chartreuse, chocolate, coral, cornflowerblue,
                cornsilk, crimson, cyan, darkblue, darkcyan,
                darkgoldenrod, darkgray, darkgrey, darkgreen,
                darkkhaki, darkmagenta, darkolivegreen, darkorange,
                darkorchid, darkred, darksalmon, darkseagreen,
                darkslateblue, darkslategray, darkslategrey,
                darkturquoise, darkviolet, deeppink, deepskyblue,
                dimgray, dimgrey, dodgerblue, firebrick,
                floralwhite, forestgreen, fuchsia, gainsboro,
                ghostwhite, gold, goldenrod, gray, grey, green,
                greenyellow, honeydew, hotpink, indianred, indigo,
                ivory, khaki, lavender, lavenderblush, lawngreen,
                lemonchiffon, lightblue, lightcoral, lightcyan,
                lightgoldenrodyellow, lightgray, lightgrey,
                lightgreen, lightpink, lightsalmon, lightseagreen,
                lightskyblue, lightslategray, lightslategrey,
                lightsteelblue, lightyellow, lime, limegreen,
                linen, magenta, maroon, mediumaquamarine,
                mediumblue, mediumorchid, mediumpurple,
                mediumseagreen, mediumslateblue, mediumspringgreen,
                mediumturquoise, mediumvioletred, midnightblue,
                mintcream, mistyrose, moccasin, navajowhite, navy,
                oldlace, olive, olivedrab, orange, orangered,
                orchid, palegoldenrod, palegreen, paleturquoise,
                palevioletred, papayawhip, peachpuff, peru, pink,
                plum, powderblue, purple, red, rosybrown,
                royalblue, rebeccapurple, saddlebrown, salmon,
                sandybrown, seagreen, seashell, sienna, silver,
                skyblue, slateblue, slategray, slategrey, snow,
                springgreen, steelblue, tan, teal, thistle, tomato,
                turquoise, violet, wheat, white, whitesmoke,
                yellow, yellowgreen
          - A list or array of any of the above

        Returns
        -------
        str|numpy.ndarray
        """
        return self["bgcolor"]

    @bgcolor.setter
    def bgcolor(self, val):
        self["bgcolor"] = val

    # bgcolorsrc
    # ----------
    @property
    def bgcolorsrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `bgcolor`.

        The 'bgcolorsrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["bgcolorsrc"]

    @bgcolorsrc.setter
    def bgcolorsrc(self, val):
        self["bgcolorsrc"] = val

    # fgcolor
    # -------
    @property
    def fgcolor(self):
        """
        When there is no colorscale sets the color of foreground
        pattern fill. Defaults to a `marker.color` background when
        `fillmode` is "replace". Otherwise, defaults to dark grey or
        white to increase contrast with the `bgcolor`.

        The 'fgcolor' property is a color and may be specified as:
          - A hex string (e.g. '#ff0000')
          - An rgb/rgba string (e.g. 'rgb(255,0,0)')
          - An hsl/hsla string (e.g. 'hsl(0,100%,50%)')
          - An hsv/hsva string (e.g. 'hsv(0,100%,100%)')
          - A named CSS color:
                aliceblue, antiquewhite, aqua, aquamarine, azure,
                beige, bisque, black, blanchedalmond, blue,
                blueviolet, brown, burlywood, cadetblue,
                chartreuse, chocolate, coral, cornflowerblue,
                cornsilk, crimson, cyan, darkblue, darkcyan,
                darkgoldenrod, darkgray, darkgrey, darkgreen,
                darkkhaki, darkmagenta, darkolivegreen, darkorange,
                darkorchid, darkred, darksalmon, darkseagreen,
                darkslateblue, darkslategray, darkslategrey,
                darkturquoise, darkviolet, deeppink, deepskyblue,
                dimgray, dimgrey, dodgerblue, firebrick,
                floralwhite, forestgreen, fuchsia, gainsboro,
                ghostwhite, gold, goldenrod, gray, grey, green,
                greenyellow, honeydew, hotpink, indianred, indigo,
                ivory, khaki, lavender, lavenderblush, lawngreen,
                lemonchiffon, lightblue, lightcoral, lightcyan,
                lightgoldenrodyellow, lightgray, lightgrey,
                lightgreen, lightpink, lightsalmon, lightseagreen,
                lightskyblue, lightslategray, lightslategrey,
                lightsteelblue, lightyellow, lime, limegreen,
                linen, magenta, maroon, mediumaquamarine,
                mediumblue, mediumorchid, mediumpurple,
                mediumseagreen, mediumslateblue, mediumspringgreen,
                mediumturquoise, mediumvioletred, midnightblue,
                mintcream, mistyrose, moccasin, navajowhite, navy,
                oldlace, olive, olivedrab, orange, orangered,
                orchid, palegoldenrod, palegreen, paleturquoise,
                palevioletred, papayawhip, peachpuff, peru, pink,
                plum, powderblue, purple, red, rosybrown,
                royalblue, rebeccapurple, saddlebrown, salmon,
                sandybrown, seagreen, seashell, sienna, silver,
                skyblue, slateblue, slategray, slategrey, snow,
                springgreen, steelblue, tan, teal, thistle, tomato,
                turquoise, violet, wheat, white, whitesmoke,
                yellow, yellowgreen
          - A list or array of any of the above

        Returns
        -------
        str|numpy.ndarray
        """
        return self["fgcolor"]

    @fgcolor.setter
    def fgcolor(self, val):
        self["fgcolor"] = val

    # fgcolorsrc
    # ----------
    @property
    def fgcolorsrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `fgcolor`.

        The 'fgcolorsrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["fgcolorsrc"]

    @fgcolorsrc.setter
    def fgcolorsrc(self, val):
        self["fgcolorsrc"] = val

    # fgopacity
    # ---------
    @property
    def fgopacity(self):
        """
        Sets the opacity of the foreground pattern fill. Defaults to a
        0.5 when `fillmode` is "overlay". Otherwise, defaults to 1.

        The 'fgopacity' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["fgopacity"]

    @fgopacity.setter
    def fgopacity(self, val):
        self["fgopacity"] = val

    # fillmode
    # --------
    @property
    def fillmode(self):
        """
        Determines whether `marker.color` should be used as a default
        to `bgcolor` or a `fgcolor`.

        The 'fillmode' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['replace', 'overlay']

        Returns
        -------
        Any
        """
        return self["fillmode"]

    @fillmode.setter
    def fillmode(self, val):
        self["fillmode"] = val

    # shape
    # -----
    @property
    def shape(self):
        """
        Sets the shape of the pattern fill. By default, no pattern is
        used for filling the area.

        The 'shape' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['', '/', '\\', 'x', '-', '|', '+', '.']
          - A tuple, list, or one-dimensional numpy array of the above

        Returns
        -------
        Any|numpy.ndarray
        """
        return self["shape"]

    @shape.setter
    def shape(self, val):
        self["shape"] = val

    # shapesrc
    # --------
    @property
    def shapesrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `shape`.

        The 'shapesrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["shapesrc"]

    @shapesrc.setter
    def shapesrc(self, val):
        self["shapesrc"] = val

    # size
    # ----
    @property
    def size(self):
        """
        Sets the size of unit squares of the pattern fill in pixels,
        which corresponds to the interval of repetition of the pattern.

        The 'size' property is a number and may be specified as:
          - An int or float in the interval [0, inf]
          - A tuple, list, or one-dimensional numpy array of the above

        Returns
        -------
        int|float|numpy.ndarray
        """
        return self["size"]

    @size.setter
    def size(self, val):
        self["size"] = val

    # sizesrc
    # -------
    @property
    def sizesrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `size`.

        The 'sizesrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["sizesrc"]

    @sizesrc.setter
    def sizesrc(self, val):
        self["sizesrc"] = val

    # solidity
    # --------
    @property
    def solidity(self):
        """
        Sets the solidity of the pattern fill. Solidity is roughly the
        fraction of the area filled by the pattern. Solidity of 0 shows
        only the background color without pattern and solidty of 1
        shows only the foreground color without pattern.

        The 'solidity' property is a number and may be specified as:
          - An int or float in the interval [0, 1]
          - A tuple, list, or one-dimensional numpy array of the above

        Returns
        -------
        int|float|numpy.ndarray
        """
        return self["solidity"]

    @solidity.setter
    def solidity(self, val):
        self["solidity"] = val

    # soliditysrc
    # -----------
    @property
    def soliditysrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `solidity`.

        The 'soliditysrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["soliditysrc"]

    @soliditysrc.setter
    def soliditysrc(self, val):
        self["soliditysrc"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        bgcolor
            When there is no colorscale sets the color of
            background pattern fill. Defaults to a `marker.color`
            background when `fillmode` is "overlay". Otherwise,
            defaults to a transparent background.
        bgcolorsrc
            Sets the source reference on Chart Studio Cloud for
            `bgcolor`.
        fgcolor
            When there is no colorscale sets the color of
            foreground pattern fill. Defaults to a `marker.color`
            background when `fillmode` is "replace". Otherwise,
            defaults to dark grey or white to increase contrast
            with the `bgcolor`.
        fgcolorsrc
            Sets the source reference on Chart Studio Cloud for
            `fgcolor`.
        fgopacity
            Sets the opacity of the foreground pattern fill.
            Defaults to a 0.5 when `fillmode` is "overlay".
            Otherwise, defaults to 1.
        fillmode
            Determines whether `marker.color` should be used as a
            default to `bgcolor` or a `fgcolor`.
        shape
            Sets the shape of the pattern fill. By default, no
            pattern is used for filling the area.
        shapesrc
            Sets the source reference on Chart Studio Cloud for
            `shape`.
        size
            Sets the size of unit squares of the pattern fill in
            pixels, which corresponds to the interval of repetition
            of the pattern.
        sizesrc
            Sets the source reference on Chart Studio Cloud for
            `size`.
        solidity
            Sets the solidity of the pattern fill. Solidity is
            roughly the fraction of the area filled by the pattern.
            Solidity of 0 shows only the background color without
            pattern and solidty of 1 shows only the foreground
            color without pattern.
        soliditysrc
            Sets the source reference on Chart Studio Cloud for
            `solidity`.
        """

    def __init__(
        self,
        arg=None,
        bgcolor=None,
        bgcolorsrc=None,
        fgcolor=None,
        fgcolorsrc=None,
        fgopacity=None,
        fillmode=None,
        shape=None,
        shapesrc=None,
        size=None,
        sizesrc=None,
        solidity=None,
        soliditysrc=None,
        **kwargs,
    ):
        """
        Construct a new Pattern object

        Sets the pattern within the marker.

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.sunburst.marker.Pattern`
        bgcolor
            When there is no colorscale sets the color of
            background pattern fill. Defaults to a `marker.color`
            background when `fillmode` is "overlay". Otherwise,
            defaults to a transparent background.
        bgcolorsrc
            Sets the source reference on Chart Studio Cloud for
            `bgcolor`.
        fgcolor
            When there is no colorscale sets the color of
            foreground pattern fill. Defaults to a `marker.color`
            background when `fillmode` is "replace". Otherwise,
            defaults to dark grey or white to increase contrast
            with the `bgcolor`.
        fgcolorsrc
            Sets the source reference on Chart Studio Cloud for
            `fgcolor`.
        fgopacity
            Sets the opacity of the foreground pattern fill.
            Defaults to a 0.5 when `fillmode` is "overlay".
            Otherwise, defaults to 1.
        fillmode
            Determines whether `marker.color` should be used as a
            default to `bgcolor` or a `fgcolor`.
        shape
            Sets the shape of the pattern fill. By default, no
            pattern is used for filling the area.
        shapesrc
            Sets the source reference on Chart Studio Cloud for
            `shape`.
        size
            Sets the size of unit squares of the pattern fill in
            pixels, which corresponds to the interval of repetition
            of the pattern.
        sizesrc
            Sets the source reference on Chart Studio Cloud for
            `size`.
        solidity
            Sets the solidity of the pattern fill. Solidity is
            roughly the fraction of the area filled by the pattern.
            Solidity of 0 shows only the background color without
            pattern and solidty of 1 shows only the foreground
            color without pattern.
        soliditysrc
            Sets the source reference on Chart Studio Cloud for
            `solidity`.

        Returns
        -------
        Pattern
        """
        super(Pattern, self).__init__("pattern")

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
The first argument to the plotly.graph_objs.sunburst.marker.Pattern
constructor must be a dict or
an instance of :class:`plotly.graph_objs.sunburst.marker.Pattern`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("bgcolor", None)
        _v = bgcolor if bgcolor is not None else _v
        if _v is not None:
            self["bgcolor"] = _v
        _v = arg.pop("bgcolorsrc", None)
        _v = bgcolorsrc if bgcolorsrc is not None else _v
        if _v is not None:
            self["bgcolorsrc"] = _v
        _v = arg.pop("fgcolor", None)
        _v = fgcolor if fgcolor is not None else _v
        if _v is not None:
            self["fgcolor"] = _v
        _v = arg.pop("fgcolorsrc", None)
        _v = fgcolorsrc if fgcolorsrc is not None else _v
        if _v is not None:
            self["fgcolorsrc"] = _v
        _v = arg.pop("fgopacity", None)
        _v = fgopacity if fgopacity is not None else _v
        if _v is not None:
            self["fgopacity"] = _v
        _v = arg.pop("fillmode", None)
        _v = fillmode if fillmode is not None else _v
        if _v is not None:
            self["fillmode"] = _v
        _v = arg.pop("shape", None)
        _v = shape if shape is not None else _v
        if _v is not None:
            self["shape"] = _v
        _v = arg.pop("shapesrc", None)
        _v = shapesrc if shapesrc is not None else _v
        if _v is not None:
            self["shapesrc"] = _v
        _v = arg.pop("size", None)
        _v = size if size is not None else _v
        if _v is not None:
            self["size"] = _v
        _v = arg.pop("sizesrc", None)
        _v = sizesrc if sizesrc is not None else _v
        if _v is not None:
            self["sizesrc"] = _v
        _v = arg.pop("solidity", None)
        _v = solidity if solidity is not None else _v
        if _v is not None:
            self["solidity"] = _v
        _v = arg.pop("soliditysrc", None)
        _v = soliditysrc if soliditysrc is not None else _v
        if _v is not None:
            self["soliditysrc"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
