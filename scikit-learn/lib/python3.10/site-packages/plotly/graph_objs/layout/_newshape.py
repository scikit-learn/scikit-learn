#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Newshape(_BaseLayoutHierarchyType):
    _parent_path_str = "layout"
    _path_str = "layout.newshape"
    _valid_props = {
        "drawdirection",
        "fillcolor",
        "fillrule",
        "label",
        "layer",
        "legend",
        "legendgroup",
        "legendgrouptitle",
        "legendrank",
        "legendwidth",
        "line",
        "name",
        "opacity",
        "showlegend",
        "visible",
    }

    @property
    def drawdirection(self):
        """
        When `dragmode` is set to "drawrect", "drawline" or
        "drawcircle" this limits the drag to be horizontal, vertical or
        diagonal. Using "diagonal" there is no limit e.g. in drawing
        lines in any direction. "ortho" limits the draw to be either
        horizontal or vertical. "horizontal" allows horizontal extend.
        "vertical" allows vertical extend.

        The 'drawdirection' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['ortho', 'horizontal', 'vertical', 'diagonal']

        Returns
        -------
        Any
        """
        return self["drawdirection"]

    @drawdirection.setter
    def drawdirection(self, val):
        self["drawdirection"] = val

    @property
    def fillcolor(self):
        """
        Sets the color filling new shapes' interior. Please note that
        if using a fillcolor with alpha greater than half, drag inside
        the active shape starts moving the shape underneath, otherwise
        a new shape could be started over.

        The 'fillcolor' property is a color and may be specified as:
          - A hex string (e.g. '#ff0000')
          - An rgb/rgba string (e.g. 'rgb(255,0,0)')
          - An hsl/hsla string (e.g. 'hsl(0,100%,50%)')
          - An hsv/hsva string (e.g. 'hsv(0,100%,100%)')
          - A named CSS color: see https://plotly.com/python/css-colors/ for a list

        Returns
        -------
        str
        """
        return self["fillcolor"]

    @fillcolor.setter
    def fillcolor(self, val):
        self["fillcolor"] = val

    @property
    def fillrule(self):
        """
        Determines the path's interior. For more info please visit
        https://developer.mozilla.org/en-
        US/docs/Web/SVG/Attribute/fill-rule

        The 'fillrule' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['evenodd', 'nonzero']

        Returns
        -------
        Any
        """
        return self["fillrule"]

    @fillrule.setter
    def fillrule(self, val):
        self["fillrule"] = val

    @property
    def label(self):
        """
        The 'label' property is an instance of Label
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.newshape.Label`
          - A dict of string/value properties that will be passed
            to the Label constructor

        Returns
        -------
        plotly.graph_objs.layout.newshape.Label
        """
        return self["label"]

    @label.setter
    def label(self, val):
        self["label"] = val

    @property
    def layer(self):
        """
        Specifies whether new shapes are drawn below gridlines
        ("below"), between gridlines and traces ("between") or above
        traces ("above").

        The 'layer' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['below', 'above', 'between']

        Returns
        -------
        Any
        """
        return self["layer"]

    @layer.setter
    def layer(self, val):
        self["layer"] = val

    @property
    def legend(self):
        """
        Sets the reference to a legend to show new shape in. References
        to these legends are "legend", "legend2", "legend3", etc.
        Settings for these legends are set in the layout, under
        `layout.legend`, `layout.legend2`, etc.

        The 'legend' property is an identifier of a particular
        subplot, of type 'legend', that may be specified as the string 'legend'
        optionally followed by an integer >= 1
        (e.g. 'legend', 'legend1', 'legend2', 'legend3', etc.)

        Returns
        -------
        str
        """
        return self["legend"]

    @legend.setter
    def legend(self, val):
        self["legend"] = val

    @property
    def legendgroup(self):
        """
        Sets the legend group for new shape. Traces and shapes part of
        the same legend group hide/show at the same time when toggling
        legend items.

        The 'legendgroup' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["legendgroup"]

    @legendgroup.setter
    def legendgroup(self, val):
        self["legendgroup"] = val

    @property
    def legendgrouptitle(self):
        """
        The 'legendgrouptitle' property is an instance of Legendgrouptitle
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.newshape.Legendgrouptitle`
          - A dict of string/value properties that will be passed
            to the Legendgrouptitle constructor

        Returns
        -------
        plotly.graph_objs.layout.newshape.Legendgrouptitle
        """
        return self["legendgrouptitle"]

    @legendgrouptitle.setter
    def legendgrouptitle(self, val):
        self["legendgrouptitle"] = val

    @property
    def legendrank(self):
        """
        Sets the legend rank for new shape. Items and groups with
        smaller ranks are presented on top/left side while with
        "reversed" `legend.traceorder` they are on bottom/right side.
        The default legendrank is 1000, so that you can use ranks less
        than 1000 to place certain items before all unranked items, and
        ranks greater than 1000 to go after all unranked items.

        The 'legendrank' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["legendrank"]

    @legendrank.setter
    def legendrank(self, val):
        self["legendrank"] = val

    @property
    def legendwidth(self):
        """
        Sets the width (in px or fraction) of the legend for new shape.

        The 'legendwidth' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["legendwidth"]

    @legendwidth.setter
    def legendwidth(self, val):
        self["legendwidth"] = val

    @property
    def line(self):
        """
        The 'line' property is an instance of Line
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.newshape.Line`
          - A dict of string/value properties that will be passed
            to the Line constructor

        Returns
        -------
        plotly.graph_objs.layout.newshape.Line
        """
        return self["line"]

    @line.setter
    def line(self, val):
        self["line"] = val

    @property
    def name(self):
        """
        Sets new shape name. The name appears as the legend item.

        The 'name' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["name"]

    @name.setter
    def name(self, val):
        self["name"] = val

    @property
    def opacity(self):
        """
        Sets the opacity of new shapes.

        The 'opacity' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["opacity"]

    @opacity.setter
    def opacity(self, val):
        self["opacity"] = val

    @property
    def showlegend(self):
        """
        Determines whether or not new shape is shown in the legend.

        The 'showlegend' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["showlegend"]

    @showlegend.setter
    def showlegend(self, val):
        self["showlegend"] = val

    @property
    def visible(self):
        """
        Determines whether or not new shape is visible. If
        "legendonly", the shape is not drawn, but can appear as a
        legend item (provided that the legend itself is visible).

        The 'visible' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                [True, False, 'legendonly']

        Returns
        -------
        Any
        """
        return self["visible"]

    @visible.setter
    def visible(self, val):
        self["visible"] = val

    @property
    def _prop_descriptions(self):
        return """\
        drawdirection
            When `dragmode` is set to "drawrect", "drawline" or
            "drawcircle" this limits the drag to be horizontal,
            vertical or diagonal. Using "diagonal" there is no
            limit e.g. in drawing lines in any direction. "ortho"
            limits the draw to be either horizontal or vertical.
            "horizontal" allows horizontal extend. "vertical"
            allows vertical extend.
        fillcolor
            Sets the color filling new shapes' interior. Please
            note that if using a fillcolor with alpha greater than
            half, drag inside the active shape starts moving the
            shape underneath, otherwise a new shape could be
            started over.
        fillrule
            Determines the path's interior. For more info please
            visit https://developer.mozilla.org/en-
            US/docs/Web/SVG/Attribute/fill-rule
        label
            :class:`plotly.graph_objects.layout.newshape.Label`
            instance or dict with compatible properties
        layer
            Specifies whether new shapes are drawn below gridlines
            ("below"), between gridlines and traces ("between") or
            above traces ("above").
        legend
            Sets the reference to a legend to show new shape in.
            References to these legends are "legend", "legend2",
            "legend3", etc. Settings for these legends are set in
            the layout, under `layout.legend`, `layout.legend2`,
            etc.
        legendgroup
            Sets the legend group for new shape. Traces and shapes
            part of the same legend group hide/show at the same
            time when toggling legend items.
        legendgrouptitle
            :class:`plotly.graph_objects.layout.newshape.Legendgrou
            ptitle` instance or dict with compatible properties
        legendrank
            Sets the legend rank for new shape. Items and groups
            with smaller ranks are presented on top/left side while
            with "reversed" `legend.traceorder` they are on
            bottom/right side. The default legendrank is 1000, so
            that you can use ranks less than 1000 to place certain
            items before all unranked items, and ranks greater than
            1000 to go after all unranked items.
        legendwidth
            Sets the width (in px or fraction) of the legend for
            new shape.
        line
            :class:`plotly.graph_objects.layout.newshape.Line`
            instance or dict with compatible properties
        name
            Sets new shape name. The name appears as the legend
            item.
        opacity
            Sets the opacity of new shapes.
        showlegend
            Determines whether or not new shape is shown in the
            legend.
        visible
            Determines whether or not new shape is visible. If
            "legendonly", the shape is not drawn, but can appear as
            a legend item (provided that the legend itself is
            visible).
        """

    def __init__(
        self,
        arg=None,
        drawdirection=None,
        fillcolor=None,
        fillrule=None,
        label=None,
        layer=None,
        legend=None,
        legendgroup=None,
        legendgrouptitle=None,
        legendrank=None,
        legendwidth=None,
        line=None,
        name=None,
        opacity=None,
        showlegend=None,
        visible=None,
        **kwargs,
    ):
        """
        Construct a new Newshape object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.layout.Newshape`
        drawdirection
            When `dragmode` is set to "drawrect", "drawline" or
            "drawcircle" this limits the drag to be horizontal,
            vertical or diagonal. Using "diagonal" there is no
            limit e.g. in drawing lines in any direction. "ortho"
            limits the draw to be either horizontal or vertical.
            "horizontal" allows horizontal extend. "vertical"
            allows vertical extend.
        fillcolor
            Sets the color filling new shapes' interior. Please
            note that if using a fillcolor with alpha greater than
            half, drag inside the active shape starts moving the
            shape underneath, otherwise a new shape could be
            started over.
        fillrule
            Determines the path's interior. For more info please
            visit https://developer.mozilla.org/en-
            US/docs/Web/SVG/Attribute/fill-rule
        label
            :class:`plotly.graph_objects.layout.newshape.Label`
            instance or dict with compatible properties
        layer
            Specifies whether new shapes are drawn below gridlines
            ("below"), between gridlines and traces ("between") or
            above traces ("above").
        legend
            Sets the reference to a legend to show new shape in.
            References to these legends are "legend", "legend2",
            "legend3", etc. Settings for these legends are set in
            the layout, under `layout.legend`, `layout.legend2`,
            etc.
        legendgroup
            Sets the legend group for new shape. Traces and shapes
            part of the same legend group hide/show at the same
            time when toggling legend items.
        legendgrouptitle
            :class:`plotly.graph_objects.layout.newshape.Legendgrou
            ptitle` instance or dict with compatible properties
        legendrank
            Sets the legend rank for new shape. Items and groups
            with smaller ranks are presented on top/left side while
            with "reversed" `legend.traceorder` they are on
            bottom/right side. The default legendrank is 1000, so
            that you can use ranks less than 1000 to place certain
            items before all unranked items, and ranks greater than
            1000 to go after all unranked items.
        legendwidth
            Sets the width (in px or fraction) of the legend for
            new shape.
        line
            :class:`plotly.graph_objects.layout.newshape.Line`
            instance or dict with compatible properties
        name
            Sets new shape name. The name appears as the legend
            item.
        opacity
            Sets the opacity of new shapes.
        showlegend
            Determines whether or not new shape is shown in the
            legend.
        visible
            Determines whether or not new shape is visible. If
            "legendonly", the shape is not drawn, but can appear as
            a legend item (provided that the legend itself is
            visible).

        Returns
        -------
        Newshape
        """
        super().__init__("newshape")
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
The first argument to the plotly.graph_objs.layout.Newshape
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.Newshape`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("drawdirection", arg, drawdirection)
        self._set_property("fillcolor", arg, fillcolor)
        self._set_property("fillrule", arg, fillrule)
        self._set_property("label", arg, label)
        self._set_property("layer", arg, layer)
        self._set_property("legend", arg, legend)
        self._set_property("legendgroup", arg, legendgroup)
        self._set_property("legendgrouptitle", arg, legendgrouptitle)
        self._set_property("legendrank", arg, legendrank)
        self._set_property("legendwidth", arg, legendwidth)
        self._set_property("line", arg, line)
        self._set_property("name", arg, name)
        self._set_property("opacity", arg, opacity)
        self._set_property("showlegend", arg, showlegend)
        self._set_property("visible", arg, visible)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
