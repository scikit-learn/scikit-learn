#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Tickfont(_BaseLayoutHierarchyType):
    _parent_path_str = "layout.ternary.baxis"
    _path_str = "layout.ternary.baxis.tickfont"
    _valid_props = {
        "color",
        "family",
        "lineposition",
        "shadow",
        "size",
        "style",
        "textcase",
        "variant",
        "weight",
    }

    @property
    def color(self):
        """
        The 'color' property is a color and may be specified as:
          - A hex string (e.g. '#ff0000')
          - An rgb/rgba string (e.g. 'rgb(255,0,0)')
          - An hsl/hsla string (e.g. 'hsl(0,100%,50%)')
          - An hsv/hsva string (e.g. 'hsv(0,100%,100%)')
          - A named CSS color: see https://plotly.com/python/css-colors/ for a list

        Returns
        -------
        str
        """
        return self["color"]

    @color.setter
    def color(self, val):
        self["color"] = val

    @property
    def family(self):
        """
        HTML font family - the typeface that will be applied by the web
        browser. The web browser can only apply a font if it is
        available on the system where it runs. Provide multiple font
        families, separated by commas, to indicate the order in which
        to apply fonts if they aren't available.

        The 'family' property is a string and must be specified as:
          - A non-empty string

        Returns
        -------
        str
        """
        return self["family"]

    @family.setter
    def family(self, val):
        self["family"] = val

    @property
    def lineposition(self):
        """
        Sets the kind of decoration line(s) with text, such as an
        "under", "over" or "through" as well as combinations e.g.
        "under+over", etc.

        The 'lineposition' property is a flaglist and may be specified
        as a string containing:
          - Any combination of ['under', 'over', 'through'] joined with '+' characters
            (e.g. 'under+over')
            OR exactly one of ['none'] (e.g. 'none')

        Returns
        -------
        Any
        """
        return self["lineposition"]

    @lineposition.setter
    def lineposition(self, val):
        self["lineposition"] = val

    @property
    def shadow(self):
        """
        Sets the shape and color of the shadow behind text. "auto"
        places minimal shadow and applies contrast text font color. See
        https://developer.mozilla.org/en-US/docs/Web/CSS/text-shadow
        for additional options.

        The 'shadow' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["shadow"]

    @shadow.setter
    def shadow(self, val):
        self["shadow"] = val

    @property
    def size(self):
        """
        The 'size' property is a number and may be specified as:
          - An int or float in the interval [1, inf]

        Returns
        -------
        int|float
        """
        return self["size"]

    @size.setter
    def size(self, val):
        self["size"] = val

    @property
    def style(self):
        """
        Sets whether a font should be styled with a normal or italic
        face from its family.

        The 'style' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['normal', 'italic']

        Returns
        -------
        Any
        """
        return self["style"]

    @style.setter
    def style(self, val):
        self["style"] = val

    @property
    def textcase(self):
        """
        Sets capitalization of text. It can be used to make text appear
        in all-uppercase or all-lowercase, or with each word
        capitalized.

        The 'textcase' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['normal', 'word caps', 'upper', 'lower']

        Returns
        -------
        Any
        """
        return self["textcase"]

    @textcase.setter
    def textcase(self, val):
        self["textcase"] = val

    @property
    def variant(self):
        """
        Sets the variant of the font.

        The 'variant' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['normal', 'small-caps', 'all-small-caps',
                'all-petite-caps', 'petite-caps', 'unicase']

        Returns
        -------
        Any
        """
        return self["variant"]

    @variant.setter
    def variant(self, val):
        self["variant"] = val

    @property
    def weight(self):
        """
        Sets the weight (or boldness) of the font.

        The 'weight' property is a integer and may be specified as:
          - An int (or float that will be cast to an int)
            in the interval [1, 1000]
            OR exactly one of ['normal', 'bold'] (e.g. 'bold')

        Returns
        -------
        int
        """
        return self["weight"]

    @weight.setter
    def weight(self, val):
        self["weight"] = val

    @property
    def _prop_descriptions(self):
        return """\
        color

        family
            HTML font family - the typeface that will be applied by
            the web browser. The web browser can only apply a font
            if it is available on the system where it runs. Provide
            multiple font families, separated by commas, to
            indicate the order in which to apply fonts if they
            aren't available.
        lineposition
            Sets the kind of decoration line(s) with text, such as
            an "under", "over" or "through" as well as combinations
            e.g. "under+over", etc.
        shadow
            Sets the shape and color of the shadow behind text.
            "auto" places minimal shadow and applies contrast text
            font color. See https://developer.mozilla.org/en-
            US/docs/Web/CSS/text-shadow for additional options.
        size

        style
            Sets whether a font should be styled with a normal or
            italic face from its family.
        textcase
            Sets capitalization of text. It can be used to make
            text appear in all-uppercase or all-lowercase, or with
            each word capitalized.
        variant
            Sets the variant of the font.
        weight
            Sets the weight (or boldness) of the font.
        """

    def __init__(
        self,
        arg=None,
        color=None,
        family=None,
        lineposition=None,
        shadow=None,
        size=None,
        style=None,
        textcase=None,
        variant=None,
        weight=None,
        **kwargs,
    ):
        """
        Construct a new Tickfont object

        Sets the tick font.

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of :class:`plotly.graph_objs.layout.ternary
            .baxis.Tickfont`
        color

        family
            HTML font family - the typeface that will be applied by
            the web browser. The web browser can only apply a font
            if it is available on the system where it runs. Provide
            multiple font families, separated by commas, to
            indicate the order in which to apply fonts if they
            aren't available.
        lineposition
            Sets the kind of decoration line(s) with text, such as
            an "under", "over" or "through" as well as combinations
            e.g. "under+over", etc.
        shadow
            Sets the shape and color of the shadow behind text.
            "auto" places minimal shadow and applies contrast text
            font color. See https://developer.mozilla.org/en-
            US/docs/Web/CSS/text-shadow for additional options.
        size

        style
            Sets whether a font should be styled with a normal or
            italic face from its family.
        textcase
            Sets capitalization of text. It can be used to make
            text appear in all-uppercase or all-lowercase, or with
            each word capitalized.
        variant
            Sets the variant of the font.
        weight
            Sets the weight (or boldness) of the font.

        Returns
        -------
        Tickfont
        """
        super().__init__("tickfont")
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
The first argument to the plotly.graph_objs.layout.ternary.baxis.Tickfont
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.ternary.baxis.Tickfont`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("color", arg, color)
        self._set_property("family", arg, family)
        self._set_property("lineposition", arg, lineposition)
        self._set_property("shadow", arg, shadow)
        self._set_property("size", arg, size)
        self._set_property("style", arg, style)
        self._set_property("textcase", arg, textcase)
        self._set_property("variant", arg, variant)
        self._set_property("weight", arg, weight)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
