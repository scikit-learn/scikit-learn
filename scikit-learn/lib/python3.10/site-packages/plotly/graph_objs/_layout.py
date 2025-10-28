#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseLayoutType as _BaseLayoutType
import copy as _copy


class Layout(_BaseLayoutType):
    _subplotid_prop_names = [
        "coloraxis",
        "geo",
        "legend",
        "map",
        "mapbox",
        "polar",
        "scene",
        "smith",
        "ternary",
        "xaxis",
        "yaxis",
    ]

    import re

    _subplotid_prop_re = re.compile("^(" + "|".join(_subplotid_prop_names) + r")(\d+)$")

    @property
    def _subplotid_validators(self):
        """
        dict of validator classes for each subplot type

        Returns
        -------
        dict
        """
        from plotly.validator_cache import ValidatorCache

        return {
            "coloraxis": ValidatorCache.get_validator("layout", "coloraxis"),
            "geo": ValidatorCache.get_validator("layout", "geo"),
            "legend": ValidatorCache.get_validator("layout", "legend"),
            "map": ValidatorCache.get_validator("layout", "map"),
            "mapbox": ValidatorCache.get_validator("layout", "mapbox"),
            "polar": ValidatorCache.get_validator("layout", "polar"),
            "scene": ValidatorCache.get_validator("layout", "scene"),
            "smith": ValidatorCache.get_validator("layout", "smith"),
            "ternary": ValidatorCache.get_validator("layout", "ternary"),
            "xaxis": ValidatorCache.get_validator("layout", "xaxis"),
            "yaxis": ValidatorCache.get_validator("layout", "yaxis"),
        }

    def _subplot_re_match(self, prop):
        return self._subplotid_prop_re.match(prop)

    _parent_path_str = ""
    _path_str = "layout"
    _valid_props = {
        "activeselection",
        "activeshape",
        "annotationdefaults",
        "annotations",
        "autosize",
        "autotypenumbers",
        "barcornerradius",
        "bargap",
        "bargroupgap",
        "barmode",
        "barnorm",
        "boxgap",
        "boxgroupgap",
        "boxmode",
        "calendar",
        "clickmode",
        "coloraxis",
        "colorscale",
        "colorway",
        "computed",
        "datarevision",
        "dragmode",
        "editrevision",
        "extendfunnelareacolors",
        "extendiciclecolors",
        "extendpiecolors",
        "extendsunburstcolors",
        "extendtreemapcolors",
        "font",
        "funnelareacolorway",
        "funnelgap",
        "funnelgroupgap",
        "funnelmode",
        "geo",
        "grid",
        "height",
        "hiddenlabels",
        "hiddenlabelssrc",
        "hidesources",
        "hoverdistance",
        "hoverlabel",
        "hovermode",
        "hoversubplots",
        "iciclecolorway",
        "imagedefaults",
        "images",
        "legend",
        "map",
        "mapbox",
        "margin",
        "meta",
        "metasrc",
        "minreducedheight",
        "minreducedwidth",
        "modebar",
        "newselection",
        "newshape",
        "paper_bgcolor",
        "piecolorway",
        "plot_bgcolor",
        "polar",
        "scattergap",
        "scattermode",
        "scene",
        "selectdirection",
        "selectiondefaults",
        "selectionrevision",
        "selections",
        "separators",
        "shapedefaults",
        "shapes",
        "showlegend",
        "sliderdefaults",
        "sliders",
        "smith",
        "spikedistance",
        "sunburstcolorway",
        "template",
        "ternary",
        "title",
        "transition",
        "treemapcolorway",
        "uirevision",
        "uniformtext",
        "updatemenudefaults",
        "updatemenus",
        "violingap",
        "violingroupgap",
        "violinmode",
        "waterfallgap",
        "waterfallgroupgap",
        "waterfallmode",
        "width",
        "xaxis",
        "yaxis",
    }

    @property
    def activeselection(self):
        """
        The 'activeselection' property is an instance of Activeselection
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Activeselection`
          - A dict of string/value properties that will be passed
            to the Activeselection constructor

        Returns
        -------
        plotly.graph_objs.layout.Activeselection
        """
        return self["activeselection"]

    @activeselection.setter
    def activeselection(self, val):
        self["activeselection"] = val

    @property
    def activeshape(self):
        """
        The 'activeshape' property is an instance of Activeshape
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Activeshape`
          - A dict of string/value properties that will be passed
            to the Activeshape constructor

        Returns
        -------
        plotly.graph_objs.layout.Activeshape
        """
        return self["activeshape"]

    @activeshape.setter
    def activeshape(self, val):
        self["activeshape"] = val

    @property
    def annotations(self):
        """
        The 'annotations' property is a tuple of instances of
        Annotation that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.Annotation
          - A list or tuple of dicts of string/value properties that
            will be passed to the Annotation constructor

        Returns
        -------
        tuple[plotly.graph_objs.layout.Annotation]
        """
        return self["annotations"]

    @annotations.setter
    def annotations(self, val):
        self["annotations"] = val

    @property
    def annotationdefaults(self):
        """
        When used in a template (as
        layout.template.layout.annotationdefaults), sets the default
        property values to use for elements of layout.annotations

        The 'annotationdefaults' property is an instance of Annotation
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Annotation`
          - A dict of string/value properties that will be passed
            to the Annotation constructor

        Returns
        -------
        plotly.graph_objs.layout.Annotation
        """
        return self["annotationdefaults"]

    @annotationdefaults.setter
    def annotationdefaults(self, val):
        self["annotationdefaults"] = val

    @property
    def autosize(self):
        """
        Determines whether or not a layout width or height that has
        been left undefined by the user is initialized on each
        relayout. Note that, regardless of this attribute, an undefined
        layout width or height is always initialized on the first call
        to plot.

        The 'autosize' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["autosize"]

    @autosize.setter
    def autosize(self, val):
        self["autosize"] = val

    @property
    def autotypenumbers(self):
        """
        Using "strict" a numeric string in trace data is not converted
        to a number. Using *convert types* a numeric string in trace
        data may be treated as a number during automatic axis `type`
        detection. This is the default value; however it could be
        overridden for individual axes.

        The 'autotypenumbers' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['convert types', 'strict']

        Returns
        -------
        Any
        """
        return self["autotypenumbers"]

    @autotypenumbers.setter
    def autotypenumbers(self, val):
        self["autotypenumbers"] = val

    @property
    def barcornerradius(self):
        """
        Sets the rounding of bar corners. May be an integer number of
        pixels, or a percentage of bar width (as a string ending in %).

        The 'barcornerradius' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["barcornerradius"]

    @barcornerradius.setter
    def barcornerradius(self, val):
        self["barcornerradius"] = val

    @property
    def bargap(self):
        """
        Sets the gap (in plot fraction) between bars of adjacent
        location coordinates.

        The 'bargap' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["bargap"]

    @bargap.setter
    def bargap(self, val):
        self["bargap"] = val

    @property
    def bargroupgap(self):
        """
        Sets the gap (in plot fraction) between bars of the same
        location coordinate.

        The 'bargroupgap' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["bargroupgap"]

    @bargroupgap.setter
    def bargroupgap(self, val):
        self["bargroupgap"] = val

    @property
    def barmode(self):
        """
        Determines how bars at the same location coordinate are
        displayed on the graph. With "stack", the bars are stacked on
        top of one another With "relative", the bars are stacked on top
        of one another, with negative values below the axis, positive
        values above With "group", the bars are plotted next to one
        another centered around the shared location. With "overlay",
        the bars are plotted over one another, you might need to reduce
        "opacity" to see multiple bars.

        The 'barmode' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['stack', 'group', 'overlay', 'relative']

        Returns
        -------
        Any
        """
        return self["barmode"]

    @barmode.setter
    def barmode(self, val):
        self["barmode"] = val

    @property
    def barnorm(self):
        """
        Sets the normalization for bar traces on the graph. With
        "fraction", the value of each bar is divided by the sum of all
        values at that location coordinate. "percent" is the same but
        multiplied by 100 to show percentages.

        The 'barnorm' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['', 'fraction', 'percent']

        Returns
        -------
        Any
        """
        return self["barnorm"]

    @barnorm.setter
    def barnorm(self, val):
        self["barnorm"] = val

    @property
    def boxgap(self):
        """
        Sets the gap (in plot fraction) between boxes of adjacent
        location coordinates. Has no effect on traces that have "width"
        set.

        The 'boxgap' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["boxgap"]

    @boxgap.setter
    def boxgap(self, val):
        self["boxgap"] = val

    @property
    def boxgroupgap(self):
        """
        Sets the gap (in plot fraction) between boxes of the same
        location coordinate. Has no effect on traces that have "width"
        set.

        The 'boxgroupgap' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["boxgroupgap"]

    @boxgroupgap.setter
    def boxgroupgap(self, val):
        self["boxgroupgap"] = val

    @property
    def boxmode(self):
        """
        Determines how boxes at the same location coordinate are
        displayed on the graph. If "group", the boxes are plotted next
        to one another centered around the shared location. If
        "overlay", the boxes are plotted over one another, you might
        need to set "opacity" to see them multiple boxes. Has no effect
        on traces that have "width" set.

        The 'boxmode' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['group', 'overlay']

        Returns
        -------
        Any
        """
        return self["boxmode"]

    @boxmode.setter
    def boxmode(self, val):
        self["boxmode"] = val

    @property
    def calendar(self):
        """
        Sets the default calendar system to use for interpreting and
        displaying dates throughout the plot.

        The 'calendar' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['chinese', 'coptic', 'discworld', 'ethiopian',
                'gregorian', 'hebrew', 'islamic', 'jalali', 'julian',
                'mayan', 'nanakshahi', 'nepali', 'persian', 'taiwan',
                'thai', 'ummalqura']

        Returns
        -------
        Any
        """
        return self["calendar"]

    @calendar.setter
    def calendar(self, val):
        self["calendar"] = val

    @property
    def clickmode(self):
        """
        Determines the mode of single click interactions. "event" is
        the default value and emits the `plotly_click` event. In
        addition this mode emits the `plotly_selected` event in drag
        modes "lasso" and "select", but with no event data attached
        (kept for compatibility reasons). The "select" flag enables
        selecting single data points via click. This mode also supports
        persistent selections, meaning that pressing Shift while
        clicking, adds to / subtracts from an existing selection.
        "select" with `hovermode`: "x" can be confusing, consider
        explicitly setting `hovermode`: "closest" when using this
        feature. Selection events are sent accordingly as long as
        "event" flag is set as well. When the "event" flag is missing,
        `plotly_click` and `plotly_selected` events are not fired.

        The 'clickmode' property is a flaglist and may be specified
        as a string containing:
          - Any combination of ['event', 'select'] joined with '+' characters
            (e.g. 'event+select')
            OR exactly one of ['none'] (e.g. 'none')

        Returns
        -------
        Any
        """
        return self["clickmode"]

    @clickmode.setter
    def clickmode(self, val):
        self["clickmode"] = val

    @property
    def coloraxis(self):
        """
        The 'coloraxis' property is an instance of Coloraxis
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Coloraxis`
          - A dict of string/value properties that will be passed
            to the Coloraxis constructor

        Returns
        -------
        plotly.graph_objs.layout.Coloraxis
        """
        return self["coloraxis"]

    @coloraxis.setter
    def coloraxis(self, val):
        self["coloraxis"] = val

    @property
    def colorscale(self):
        """
        The 'colorscale' property is an instance of Colorscale
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Colorscale`
          - A dict of string/value properties that will be passed
            to the Colorscale constructor

        Returns
        -------
        plotly.graph_objs.layout.Colorscale
        """
        return self["colorscale"]

    @colorscale.setter
    def colorscale(self, val):
        self["colorscale"] = val

    @property
    def colorway(self):
        """
        Sets the default trace colors.

        The 'colorway' property is a colorlist that may be specified
        as a tuple, list, one-dimensional numpy array, or pandas Series of valid
        color strings

        Returns
        -------
        list
        """
        return self["colorway"]

    @colorway.setter
    def colorway(self, val):
        self["colorway"] = val

    @property
    def computed(self):
        """
        Placeholder for exporting automargin-impacting values namely
        `margin.t`, `margin.b`, `margin.l` and `margin.r` in "full-
        json" mode.

        The 'computed' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["computed"]

    @computed.setter
    def computed(self, val):
        self["computed"] = val

    @property
    def datarevision(self):
        """
        If provided, a changed value tells `Plotly.react` that one or
        more data arrays has changed. This way you can modify arrays
        in-place rather than making a complete new copy for an
        incremental change. If NOT provided, `Plotly.react` assumes
        that data arrays are being treated as immutable, thus any data
        array with a different identity from its predecessor contains
        new data.

        The 'datarevision' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["datarevision"]

    @datarevision.setter
    def datarevision(self, val):
        self["datarevision"] = val

    @property
    def dragmode(self):
        """
        Determines the mode of drag interactions. "select" and "lasso"
        apply only to scatter traces with markers or text. "orbit" and
        "turntable" apply only to 3D scenes.

        The 'dragmode' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['zoom', 'pan', 'select', 'lasso', 'drawclosedpath',
                'drawopenpath', 'drawline', 'drawrect', 'drawcircle',
                'orbit', 'turntable', False]

        Returns
        -------
        Any
        """
        return self["dragmode"]

    @dragmode.setter
    def dragmode(self, val):
        self["dragmode"] = val

    @property
    def editrevision(self):
        """
        Controls persistence of user-driven changes in `editable: true`
        configuration, other than trace names and axis titles. Defaults
        to `layout.uirevision`.

        The 'editrevision' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["editrevision"]

    @editrevision.setter
    def editrevision(self, val):
        self["editrevision"] = val

    @property
    def extendfunnelareacolors(self):
        """
        If `true`, the funnelarea slice colors (whether given by
        `funnelareacolorway` or inherited from `colorway`) will be
        extended to three times its original length by first repeating
        every color 20% lighter then each color 20% darker. This is
        intended to reduce the likelihood of reusing the same color
        when you have many slices, but you can set `false` to disable.
        Colors provided in the trace, using `marker.colors`, are never
        extended.

        The 'extendfunnelareacolors' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["extendfunnelareacolors"]

    @extendfunnelareacolors.setter
    def extendfunnelareacolors(self, val):
        self["extendfunnelareacolors"] = val

    @property
    def extendiciclecolors(self):
        """
        If `true`, the icicle slice colors (whether given by
        `iciclecolorway` or inherited from `colorway`) will be extended
        to three times its original length by first repeating every
        color 20% lighter then each color 20% darker. This is intended
        to reduce the likelihood of reusing the same color when you
        have many slices, but you can set `false` to disable. Colors
        provided in the trace, using `marker.colors`, are never
        extended.

        The 'extendiciclecolors' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["extendiciclecolors"]

    @extendiciclecolors.setter
    def extendiciclecolors(self, val):
        self["extendiciclecolors"] = val

    @property
    def extendpiecolors(self):
        """
        If `true`, the pie slice colors (whether given by `piecolorway`
        or inherited from `colorway`) will be extended to three times
        its original length by first repeating every color 20% lighter
        then each color 20% darker. This is intended to reduce the
        likelihood of reusing the same color when you have many slices,
        but you can set `false` to disable. Colors provided in the
        trace, using `marker.colors`, are never extended.

        The 'extendpiecolors' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["extendpiecolors"]

    @extendpiecolors.setter
    def extendpiecolors(self, val):
        self["extendpiecolors"] = val

    @property
    def extendsunburstcolors(self):
        """
        If `true`, the sunburst slice colors (whether given by
        `sunburstcolorway` or inherited from `colorway`) will be
        extended to three times its original length by first repeating
        every color 20% lighter then each color 20% darker. This is
        intended to reduce the likelihood of reusing the same color
        when you have many slices, but you can set `false` to disable.
        Colors provided in the trace, using `marker.colors`, are never
        extended.

        The 'extendsunburstcolors' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["extendsunburstcolors"]

    @extendsunburstcolors.setter
    def extendsunburstcolors(self, val):
        self["extendsunburstcolors"] = val

    @property
    def extendtreemapcolors(self):
        """
        If `true`, the treemap slice colors (whether given by
        `treemapcolorway` or inherited from `colorway`) will be
        extended to three times its original length by first repeating
        every color 20% lighter then each color 20% darker. This is
        intended to reduce the likelihood of reusing the same color
        when you have many slices, but you can set `false` to disable.
        Colors provided in the trace, using `marker.colors`, are never
        extended.

        The 'extendtreemapcolors' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["extendtreemapcolors"]

    @extendtreemapcolors.setter
    def extendtreemapcolors(self, val):
        self["extendtreemapcolors"] = val

    @property
    def font(self):
        """
        Sets the global font. Note that fonts used in traces and other
        layout components inherit from the global font.

        The 'font' property is an instance of Font
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Font`
          - A dict of string/value properties that will be passed
            to the Font constructor

        Returns
        -------
        plotly.graph_objs.layout.Font
        """
        return self["font"]

    @font.setter
    def font(self, val):
        self["font"] = val

    @property
    def funnelareacolorway(self):
        """
        Sets the default funnelarea slice colors. Defaults to the main
        `colorway` used for trace colors. If you specify a new list
        here it can still be extended with lighter and darker colors,
        see `extendfunnelareacolors`.

        The 'funnelareacolorway' property is a colorlist that may be specified
        as a tuple, list, one-dimensional numpy array, or pandas Series of valid
        color strings

        Returns
        -------
        list
        """
        return self["funnelareacolorway"]

    @funnelareacolorway.setter
    def funnelareacolorway(self, val):
        self["funnelareacolorway"] = val

    @property
    def funnelgap(self):
        """
        Sets the gap (in plot fraction) between bars of adjacent
        location coordinates.

        The 'funnelgap' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["funnelgap"]

    @funnelgap.setter
    def funnelgap(self, val):
        self["funnelgap"] = val

    @property
    def funnelgroupgap(self):
        """
        Sets the gap (in plot fraction) between bars of the same
        location coordinate.

        The 'funnelgroupgap' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["funnelgroupgap"]

    @funnelgroupgap.setter
    def funnelgroupgap(self, val):
        self["funnelgroupgap"] = val

    @property
    def funnelmode(self):
        """
        Determines how bars at the same location coordinate are
        displayed on the graph. With "stack", the bars are stacked on
        top of one another With "group", the bars are plotted next to
        one another centered around the shared location. With
        "overlay", the bars are plotted over one another, you might
        need to reduce "opacity" to see multiple bars.

        The 'funnelmode' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['stack', 'group', 'overlay']

        Returns
        -------
        Any
        """
        return self["funnelmode"]

    @funnelmode.setter
    def funnelmode(self, val):
        self["funnelmode"] = val

    @property
    def geo(self):
        """
        The 'geo' property is an instance of Geo
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Geo`
          - A dict of string/value properties that will be passed
            to the Geo constructor

        Returns
        -------
        plotly.graph_objs.layout.Geo
        """
        return self["geo"]

    @geo.setter
    def geo(self, val):
        self["geo"] = val

    @property
    def grid(self):
        """
        The 'grid' property is an instance of Grid
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Grid`
          - A dict of string/value properties that will be passed
            to the Grid constructor

        Returns
        -------
        plotly.graph_objs.layout.Grid
        """
        return self["grid"]

    @grid.setter
    def grid(self, val):
        self["grid"] = val

    @property
    def height(self):
        """
        Sets the plot's height (in px).

        The 'height' property is a number and may be specified as:
          - An int or float in the interval [10, inf]

        Returns
        -------
        int|float
        """
        return self["height"]

    @height.setter
    def height(self, val):
        self["height"] = val

    @property
    def hiddenlabels(self):
        """
        hiddenlabels is the funnelarea & pie chart analog of
        visible:'legendonly' but it can contain many labels, and can
        simultaneously hide slices from several pies/funnelarea charts

        The 'hiddenlabels' property is an array that may be specified as a tuple,
        list, numpy array, or pandas Series

        Returns
        -------
        numpy.ndarray
        """
        return self["hiddenlabels"]

    @hiddenlabels.setter
    def hiddenlabels(self, val):
        self["hiddenlabels"] = val

    @property
    def hiddenlabelssrc(self):
        """
        Sets the source reference on Chart Studio Cloud for
        `hiddenlabels`.

        The 'hiddenlabelssrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["hiddenlabelssrc"]

    @hiddenlabelssrc.setter
    def hiddenlabelssrc(self, val):
        self["hiddenlabelssrc"] = val

    @property
    def hidesources(self):
        """
        Determines whether or not a text link citing the data source is
        placed at the bottom-right cored of the figure. Has only an
        effect only on graphs that have been generated via forked
        graphs from the Chart Studio Cloud (at https://chart-
        studio.plotly.com or on-premise).

        The 'hidesources' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["hidesources"]

    @hidesources.setter
    def hidesources(self, val):
        self["hidesources"] = val

    @property
    def hoverdistance(self):
        """
        Sets the default distance (in pixels) to look for data to add
        hover labels (-1 means no cutoff, 0 means no looking for data).
        This is only a real distance for hovering on point-like
        objects, like scatter points. For area-like objects (bars,
        scatter fills, etc) hovering is on inside the area and off
        outside, but these objects will not supersede hover on point-
        like objects in case of conflict.

        The 'hoverdistance' property is a integer and may be specified as:
          - An int (or float that will be cast to an int)
            in the interval [-1, 9223372036854775807]

        Returns
        -------
        int
        """
        return self["hoverdistance"]

    @hoverdistance.setter
    def hoverdistance(self, val):
        self["hoverdistance"] = val

    @property
    def hoverlabel(self):
        """
        The 'hoverlabel' property is an instance of Hoverlabel
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Hoverlabel`
          - A dict of string/value properties that will be passed
            to the Hoverlabel constructor

        Returns
        -------
        plotly.graph_objs.layout.Hoverlabel
        """
        return self["hoverlabel"]

    @hoverlabel.setter
    def hoverlabel(self, val):
        self["hoverlabel"] = val

    @property
    def hovermode(self):
        """
        Determines the mode of hover interactions. If "closest", a
        single hoverlabel will appear for the "closest" point within
        the `hoverdistance`. If "x" (or "y"), multiple hoverlabels will
        appear for multiple points at the "closest" x- (or y-)
        coordinate within the `hoverdistance`, with the caveat that no
        more than one hoverlabel will appear per trace. If *x unified*
        (or *y unified*), a single hoverlabel will appear multiple
        points at the closest x- (or y-) coordinate within the
        `hoverdistance` with the caveat that no more than one
        hoverlabel will appear per trace. In this mode, spikelines are
        enabled by default perpendicular to the specified axis. If
        false, hover interactions are disabled.

        The 'hovermode' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['x', 'y', 'closest', False, 'x unified', 'y unified']

        Returns
        -------
        Any
        """
        return self["hovermode"]

    @hovermode.setter
    def hovermode(self, val):
        self["hovermode"] = val

    @property
    def hoversubplots(self):
        """
        Determines expansion of hover effects to other subplots If
        "single" just the axis pair of the primary point is included
        without overlaying subplots. If "overlaying" all subplots using
        the main axis and occupying the same space are included. If
        "axis", also include stacked subplots using the same axis when
        `hovermode` is set to "x", *x unified*, "y" or *y unified*.

        The 'hoversubplots' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['single', 'overlaying', 'axis']

        Returns
        -------
        Any
        """
        return self["hoversubplots"]

    @hoversubplots.setter
    def hoversubplots(self, val):
        self["hoversubplots"] = val

    @property
    def iciclecolorway(self):
        """
        Sets the default icicle slice colors. Defaults to the main
        `colorway` used for trace colors. If you specify a new list
        here it can still be extended with lighter and darker colors,
        see `extendiciclecolors`.

        The 'iciclecolorway' property is a colorlist that may be specified
        as a tuple, list, one-dimensional numpy array, or pandas Series of valid
        color strings

        Returns
        -------
        list
        """
        return self["iciclecolorway"]

    @iciclecolorway.setter
    def iciclecolorway(self, val):
        self["iciclecolorway"] = val

    @property
    def images(self):
        """
        The 'images' property is a tuple of instances of
        Image that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.Image
          - A list or tuple of dicts of string/value properties that
            will be passed to the Image constructor

        Returns
        -------
        tuple[plotly.graph_objs.layout.Image]
        """
        return self["images"]

    @images.setter
    def images(self, val):
        self["images"] = val

    @property
    def imagedefaults(self):
        """
        When used in a template (as
        layout.template.layout.imagedefaults), sets the default
        property values to use for elements of layout.images

        The 'imagedefaults' property is an instance of Image
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Image`
          - A dict of string/value properties that will be passed
            to the Image constructor

        Returns
        -------
        plotly.graph_objs.layout.Image
        """
        return self["imagedefaults"]

    @imagedefaults.setter
    def imagedefaults(self, val):
        self["imagedefaults"] = val

    @property
    def legend(self):
        """
        The 'legend' property is an instance of Legend
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Legend`
          - A dict of string/value properties that will be passed
            to the Legend constructor

        Returns
        -------
        plotly.graph_objs.layout.Legend
        """
        return self["legend"]

    @legend.setter
    def legend(self, val):
        self["legend"] = val

    @property
    def map(self):
        """
        The 'map' property is an instance of Map
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Map`
          - A dict of string/value properties that will be passed
            to the Map constructor

        Returns
        -------
        plotly.graph_objs.layout.Map
        """
        return self["map"]

    @map.setter
    def map(self, val):
        self["map"] = val

    @property
    def mapbox(self):
        """
        The 'mapbox' property is an instance of Mapbox
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Mapbox`
          - A dict of string/value properties that will be passed
            to the Mapbox constructor

        Returns
        -------
        plotly.graph_objs.layout.Mapbox
        """
        return self["mapbox"]

    @mapbox.setter
    def mapbox(self, val):
        self["mapbox"] = val

    @property
    def margin(self):
        """
        The 'margin' property is an instance of Margin
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Margin`
          - A dict of string/value properties that will be passed
            to the Margin constructor

        Returns
        -------
        plotly.graph_objs.layout.Margin
        """
        return self["margin"]

    @margin.setter
    def margin(self, val):
        self["margin"] = val

    @property
    def meta(self):
        """
        Assigns extra meta information that can be used in various
        `text` attributes. Attributes such as the graph, axis and
        colorbar `title.text`, annotation `text` `trace.name` in legend
        items, `rangeselector`, `updatemenus` and `sliders` `label`
        text all support `meta`. One can access `meta` fields using
        template strings: `%{meta[i]}` where `i` is the index of the
        `meta` item in question. `meta` can also be an object for
        example `{key: value}` which can be accessed %{meta[key]}.

        The 'meta' property accepts values of any type

        Returns
        -------
        Any|numpy.ndarray
        """
        return self["meta"]

    @meta.setter
    def meta(self, val):
        self["meta"] = val

    @property
    def metasrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `meta`.

        The 'metasrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["metasrc"]

    @metasrc.setter
    def metasrc(self, val):
        self["metasrc"] = val

    @property
    def minreducedheight(self):
        """
        Minimum height of the plot with margin.automargin applied (in
        px)

        The 'minreducedheight' property is a number and may be specified as:
          - An int or float in the interval [2, inf]

        Returns
        -------
        int|float
        """
        return self["minreducedheight"]

    @minreducedheight.setter
    def minreducedheight(self, val):
        self["minreducedheight"] = val

    @property
    def minreducedwidth(self):
        """
        Minimum width of the plot with margin.automargin applied (in
        px)

        The 'minreducedwidth' property is a number and may be specified as:
          - An int or float in the interval [2, inf]

        Returns
        -------
        int|float
        """
        return self["minreducedwidth"]

    @minreducedwidth.setter
    def minreducedwidth(self, val):
        self["minreducedwidth"] = val

    @property
    def modebar(self):
        """
        The 'modebar' property is an instance of Modebar
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Modebar`
          - A dict of string/value properties that will be passed
            to the Modebar constructor

        Returns
        -------
        plotly.graph_objs.layout.Modebar
        """
        return self["modebar"]

    @modebar.setter
    def modebar(self, val):
        self["modebar"] = val

    @property
    def newselection(self):
        """
        The 'newselection' property is an instance of Newselection
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Newselection`
          - A dict of string/value properties that will be passed
            to the Newselection constructor

        Returns
        -------
        plotly.graph_objs.layout.Newselection
        """
        return self["newselection"]

    @newselection.setter
    def newselection(self, val):
        self["newselection"] = val

    @property
    def newshape(self):
        """
        The 'newshape' property is an instance of Newshape
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Newshape`
          - A dict of string/value properties that will be passed
            to the Newshape constructor

        Returns
        -------
        plotly.graph_objs.layout.Newshape
        """
        return self["newshape"]

    @newshape.setter
    def newshape(self, val):
        self["newshape"] = val

    @property
    def paper_bgcolor(self):
        """
        Sets the background color of the paper where the graph is
        drawn.

        The 'paper_bgcolor' property is a color and may be specified as:
          - A hex string (e.g. '#ff0000')
          - An rgb/rgba string (e.g. 'rgb(255,0,0)')
          - An hsl/hsla string (e.g. 'hsl(0,100%,50%)')
          - An hsv/hsva string (e.g. 'hsv(0,100%,100%)')
          - A named CSS color: see https://plotly.com/python/css-colors/ for a list

        Returns
        -------
        str
        """
        return self["paper_bgcolor"]

    @paper_bgcolor.setter
    def paper_bgcolor(self, val):
        self["paper_bgcolor"] = val

    @property
    def piecolorway(self):
        """
        Sets the default pie slice colors. Defaults to the main
        `colorway` used for trace colors. If you specify a new list
        here it can still be extended with lighter and darker colors,
        see `extendpiecolors`.

        The 'piecolorway' property is a colorlist that may be specified
        as a tuple, list, one-dimensional numpy array, or pandas Series of valid
        color strings

        Returns
        -------
        list
        """
        return self["piecolorway"]

    @piecolorway.setter
    def piecolorway(self, val):
        self["piecolorway"] = val

    @property
    def plot_bgcolor(self):
        """
        Sets the background color of the plotting area in-between x and
        y axes.

        The 'plot_bgcolor' property is a color and may be specified as:
          - A hex string (e.g. '#ff0000')
          - An rgb/rgba string (e.g. 'rgb(255,0,0)')
          - An hsl/hsla string (e.g. 'hsl(0,100%,50%)')
          - An hsv/hsva string (e.g. 'hsv(0,100%,100%)')
          - A named CSS color: see https://plotly.com/python/css-colors/ for a list

        Returns
        -------
        str
        """
        return self["plot_bgcolor"]

    @plot_bgcolor.setter
    def plot_bgcolor(self, val):
        self["plot_bgcolor"] = val

    @property
    def polar(self):
        """
        The 'polar' property is an instance of Polar
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Polar`
          - A dict of string/value properties that will be passed
            to the Polar constructor

        Returns
        -------
        plotly.graph_objs.layout.Polar
        """
        return self["polar"]

    @polar.setter
    def polar(self, val):
        self["polar"] = val

    @property
    def scattergap(self):
        """
        Sets the gap (in plot fraction) between scatter points of
        adjacent location coordinates. Defaults to `bargap`.

        The 'scattergap' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["scattergap"]

    @scattergap.setter
    def scattergap(self, val):
        self["scattergap"] = val

    @property
    def scattermode(self):
        """
        Determines how scatter points at the same location coordinate
        are displayed on the graph. With "group", the scatter points
        are plotted next to one another centered around the shared
        location. With "overlay", the scatter points are plotted over
        one another, you might need to reduce "opacity" to see multiple
        scatter points.

        The 'scattermode' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['group', 'overlay']

        Returns
        -------
        Any
        """
        return self["scattermode"]

    @scattermode.setter
    def scattermode(self, val):
        self["scattermode"] = val

    @property
    def scene(self):
        """
        The 'scene' property is an instance of Scene
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Scene`
          - A dict of string/value properties that will be passed
            to the Scene constructor

        Returns
        -------
        plotly.graph_objs.layout.Scene
        """
        return self["scene"]

    @scene.setter
    def scene(self, val):
        self["scene"] = val

    @property
    def selectdirection(self):
        """
        When `dragmode` is set to "select", this limits the selection
        of the drag to horizontal, vertical or diagonal. "h" only
        allows horizontal selection, "v" only vertical, "d" only
        diagonal and "any" sets no limit.

        The 'selectdirection' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['h', 'v', 'd', 'any']

        Returns
        -------
        Any
        """
        return self["selectdirection"]

    @selectdirection.setter
    def selectdirection(self, val):
        self["selectdirection"] = val

    @property
    def selectionrevision(self):
        """
        Controls persistence of user-driven changes in selected points
        from all traces.

        The 'selectionrevision' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["selectionrevision"]

    @selectionrevision.setter
    def selectionrevision(self, val):
        self["selectionrevision"] = val

    @property
    def selections(self):
        """
        The 'selections' property is a tuple of instances of
        Selection that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.Selection
          - A list or tuple of dicts of string/value properties that
            will be passed to the Selection constructor

        Returns
        -------
        tuple[plotly.graph_objs.layout.Selection]
        """
        return self["selections"]

    @selections.setter
    def selections(self, val):
        self["selections"] = val

    @property
    def selectiondefaults(self):
        """
        When used in a template (as
        layout.template.layout.selectiondefaults), sets the default
        property values to use for elements of layout.selections

        The 'selectiondefaults' property is an instance of Selection
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Selection`
          - A dict of string/value properties that will be passed
            to the Selection constructor

        Returns
        -------
        plotly.graph_objs.layout.Selection
        """
        return self["selectiondefaults"]

    @selectiondefaults.setter
    def selectiondefaults(self, val):
        self["selectiondefaults"] = val

    @property
    def separators(self):
        """
        Sets the decimal and thousand separators. For example, *. *
        puts a '.' before decimals and a space between thousands. In
        English locales, dflt is ".," but other locales may alter this
        default.

        The 'separators' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["separators"]

    @separators.setter
    def separators(self, val):
        self["separators"] = val

    @property
    def shapes(self):
        """
        The 'shapes' property is a tuple of instances of
        Shape that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.Shape
          - A list or tuple of dicts of string/value properties that
            will be passed to the Shape constructor

        Returns
        -------
        tuple[plotly.graph_objs.layout.Shape]
        """
        return self["shapes"]

    @shapes.setter
    def shapes(self, val):
        self["shapes"] = val

    @property
    def shapedefaults(self):
        """
        When used in a template (as
        layout.template.layout.shapedefaults), sets the default
        property values to use for elements of layout.shapes

        The 'shapedefaults' property is an instance of Shape
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Shape`
          - A dict of string/value properties that will be passed
            to the Shape constructor

        Returns
        -------
        plotly.graph_objs.layout.Shape
        """
        return self["shapedefaults"]

    @shapedefaults.setter
    def shapedefaults(self, val):
        self["shapedefaults"] = val

    @property
    def showlegend(self):
        """
        Determines whether or not a legend is drawn. Default is `true`
        if there is a trace to show and any of these: a) Two or more
        traces would by default be shown in the legend. b) One pie
        trace is shown in the legend. c) One trace is explicitly given
        with `showlegend: true`.

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
    def sliders(self):
        """
        The 'sliders' property is a tuple of instances of
        Slider that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.Slider
          - A list or tuple of dicts of string/value properties that
            will be passed to the Slider constructor

        Returns
        -------
        tuple[plotly.graph_objs.layout.Slider]
        """
        return self["sliders"]

    @sliders.setter
    def sliders(self, val):
        self["sliders"] = val

    @property
    def sliderdefaults(self):
        """
        When used in a template (as
        layout.template.layout.sliderdefaults), sets the default
        property values to use for elements of layout.sliders

        The 'sliderdefaults' property is an instance of Slider
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Slider`
          - A dict of string/value properties that will be passed
            to the Slider constructor

        Returns
        -------
        plotly.graph_objs.layout.Slider
        """
        return self["sliderdefaults"]

    @sliderdefaults.setter
    def sliderdefaults(self, val):
        self["sliderdefaults"] = val

    @property
    def smith(self):
        """
        The 'smith' property is an instance of Smith
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Smith`
          - A dict of string/value properties that will be passed
            to the Smith constructor

        Returns
        -------
        plotly.graph_objs.layout.Smith
        """
        return self["smith"]

    @smith.setter
    def smith(self, val):
        self["smith"] = val

    @property
    def spikedistance(self):
        """
        Sets the default distance (in pixels) to look for data to draw
        spikelines to (-1 means no cutoff, 0 means no looking for
        data). As with hoverdistance, distance does not apply to area-
        like objects. In addition, some objects can be hovered on but
        will not generate spikelines, such as scatter fills.

        The 'spikedistance' property is a integer and may be specified as:
          - An int (or float that will be cast to an int)
            in the interval [-1, 9223372036854775807]

        Returns
        -------
        int
        """
        return self["spikedistance"]

    @spikedistance.setter
    def spikedistance(self, val):
        self["spikedistance"] = val

    @property
    def sunburstcolorway(self):
        """
        Sets the default sunburst slice colors. Defaults to the main
        `colorway` used for trace colors. If you specify a new list
        here it can still be extended with lighter and darker colors,
        see `extendsunburstcolors`.

        The 'sunburstcolorway' property is a colorlist that may be specified
        as a tuple, list, one-dimensional numpy array, or pandas Series of valid
        color strings

        Returns
        -------
        list
        """
        return self["sunburstcolorway"]

    @sunburstcolorway.setter
    def sunburstcolorway(self, val):
        self["sunburstcolorway"] = val

    @property
    def template(self):
        """
        Default attributes to be applied to the plot. This should be a
        dict with format: `{'layout': layoutTemplate, 'data':
        {trace_type: [traceTemplate, ...], ...}}` where
        `layoutTemplate` is a dict matching the structure of
        `figure.layout` and `traceTemplate` is a dict matching the
        structure of the trace with type `trace_type` (e.g. 'scatter').
        Alternatively, this may be specified as an instance of
        plotly.graph_objs.layout.Template.  Trace templates are applied
        cyclically to traces of each type. Container arrays (eg
        `annotations`) have special handling: An object ending in
        `defaults` (eg `annotationdefaults`) is applied to each array
        item. But if an item has a `templateitemname` key we look in
        the template array for an item with matching `name` and apply
        that instead. If no matching `name` is found we mark the item
        invisible. Any named template item not referenced is appended
        to the end of the array, so this can be used to add a watermark
        annotation or a logo image, for example. To omit one of these
        items on the plot, make an item with matching
        `templateitemname` and `visible: false`.

        The 'template' property is an instance of Template
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Template`
          - A dict of string/value properties that will be passed
            to the Template constructor
          - The name of a registered template where current registered templates
            are stored in the plotly.io.templates configuration object. The names
            of all registered templates can be retrieved with:
                >>> import plotly.io as pio
                >>> list(pio.templates)  # doctest: +ELLIPSIS
                ['ggplot2', 'seaborn', 'simple_white', 'plotly', 'plotly_white', ...]

          - A string containing multiple registered template names, joined on '+'
            characters (e.g. 'template1+template2'). In this case the resulting
            template is computed by merging together the collection of registered
            templates

        Returns
        -------
        plotly.graph_objs.layout.Template
        """
        return self["template"]

    @template.setter
    def template(self, val):
        self["template"] = val

    @property
    def ternary(self):
        """
        The 'ternary' property is an instance of Ternary
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Ternary`
          - A dict of string/value properties that will be passed
            to the Ternary constructor

        Returns
        -------
        plotly.graph_objs.layout.Ternary
        """
        return self["ternary"]

    @ternary.setter
    def ternary(self, val):
        self["ternary"] = val

    @property
    def title(self):
        """
        The 'title' property is an instance of Title
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Title`
          - A dict of string/value properties that will be passed
            to the Title constructor

        Returns
        -------
        plotly.graph_objs.layout.Title
        """
        return self["title"]

    @title.setter
    def title(self, val):
        self["title"] = val

    @property
    def transition(self):
        """
        Sets transition options used during Plotly.react updates.

        The 'transition' property is an instance of Transition
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Transition`
          - A dict of string/value properties that will be passed
            to the Transition constructor

        Returns
        -------
        plotly.graph_objs.layout.Transition
        """
        return self["transition"]

    @transition.setter
    def transition(self, val):
        self["transition"] = val

    @property
    def treemapcolorway(self):
        """
        Sets the default treemap slice colors. Defaults to the main
        `colorway` used for trace colors. If you specify a new list
        here it can still be extended with lighter and darker colors,
        see `extendtreemapcolors`.

        The 'treemapcolorway' property is a colorlist that may be specified
        as a tuple, list, one-dimensional numpy array, or pandas Series of valid
        color strings

        Returns
        -------
        list
        """
        return self["treemapcolorway"]

    @treemapcolorway.setter
    def treemapcolorway(self, val):
        self["treemapcolorway"] = val

    @property
    def uirevision(self):
        """
        Used to allow user interactions with the plot to persist after
        `Plotly.react` calls that are unaware of these interactions. If
        `uirevision` is omitted, or if it is given and it changed from
        the previous `Plotly.react` call, the exact new figure is used.
        If `uirevision` is truthy and did NOT change, any attribute
        that has been affected by user interactions and did not receive
        a different value in the new figure will keep the interaction
        value. `layout.uirevision` attribute serves as the default for
        `uirevision` attributes in various sub-containers. For finer
        control you can set these sub-attributes directly. For example,
        if your app separately controls the data on the x and y axes
        you might set `xaxis.uirevision=*time*` and
        `yaxis.uirevision=*cost*`. Then if only the y data is changed,
        you can update `yaxis.uirevision=*quantity*` and the y axis
        range will reset but the x axis range will retain any user-
        driven zoom.

        The 'uirevision' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["uirevision"]

    @uirevision.setter
    def uirevision(self, val):
        self["uirevision"] = val

    @property
    def uniformtext(self):
        """
        The 'uniformtext' property is an instance of Uniformtext
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Uniformtext`
          - A dict of string/value properties that will be passed
            to the Uniformtext constructor

        Returns
        -------
        plotly.graph_objs.layout.Uniformtext
        """
        return self["uniformtext"]

    @uniformtext.setter
    def uniformtext(self, val):
        self["uniformtext"] = val

    @property
    def updatemenus(self):
        """
        The 'updatemenus' property is a tuple of instances of
        Updatemenu that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.Updatemenu
          - A list or tuple of dicts of string/value properties that
            will be passed to the Updatemenu constructor

        Returns
        -------
        tuple[plotly.graph_objs.layout.Updatemenu]
        """
        return self["updatemenus"]

    @updatemenus.setter
    def updatemenus(self, val):
        self["updatemenus"] = val

    @property
    def updatemenudefaults(self):
        """
        When used in a template (as
        layout.template.layout.updatemenudefaults), sets the default
        property values to use for elements of layout.updatemenus

        The 'updatemenudefaults' property is an instance of Updatemenu
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.Updatemenu`
          - A dict of string/value properties that will be passed
            to the Updatemenu constructor

        Returns
        -------
        plotly.graph_objs.layout.Updatemenu
        """
        return self["updatemenudefaults"]

    @updatemenudefaults.setter
    def updatemenudefaults(self, val):
        self["updatemenudefaults"] = val

    @property
    def violingap(self):
        """
        Sets the gap (in plot fraction) between violins of adjacent
        location coordinates. Has no effect on traces that have "width"
        set.

        The 'violingap' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["violingap"]

    @violingap.setter
    def violingap(self, val):
        self["violingap"] = val

    @property
    def violingroupgap(self):
        """
        Sets the gap (in plot fraction) between violins of the same
        location coordinate. Has no effect on traces that have "width"
        set.

        The 'violingroupgap' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["violingroupgap"]

    @violingroupgap.setter
    def violingroupgap(self, val):
        self["violingroupgap"] = val

    @property
    def violinmode(self):
        """
        Determines how violins at the same location coordinate are
        displayed on the graph. If "group", the violins are plotted
        next to one another centered around the shared location. If
        "overlay", the violins are plotted over one another, you might
        need to set "opacity" to see them multiple violins. Has no
        effect on traces that have "width" set.

        The 'violinmode' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['group', 'overlay']

        Returns
        -------
        Any
        """
        return self["violinmode"]

    @violinmode.setter
    def violinmode(self, val):
        self["violinmode"] = val

    @property
    def waterfallgap(self):
        """
        Sets the gap (in plot fraction) between bars of adjacent
        location coordinates.

        The 'waterfallgap' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["waterfallgap"]

    @waterfallgap.setter
    def waterfallgap(self, val):
        self["waterfallgap"] = val

    @property
    def waterfallgroupgap(self):
        """
        Sets the gap (in plot fraction) between bars of the same
        location coordinate.

        The 'waterfallgroupgap' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["waterfallgroupgap"]

    @waterfallgroupgap.setter
    def waterfallgroupgap(self, val):
        self["waterfallgroupgap"] = val

    @property
    def waterfallmode(self):
        """
        Determines how bars at the same location coordinate are
        displayed on the graph. With "group", the bars are plotted next
        to one another centered around the shared location. With
        "overlay", the bars are plotted over one another, you might
        need to reduce "opacity" to see multiple bars.

        The 'waterfallmode' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['group', 'overlay']

        Returns
        -------
        Any
        """
        return self["waterfallmode"]

    @waterfallmode.setter
    def waterfallmode(self, val):
        self["waterfallmode"] = val

    @property
    def width(self):
        """
        Sets the plot's width (in px).

        The 'width' property is a number and may be specified as:
          - An int or float in the interval [10, inf]

        Returns
        -------
        int|float
        """
        return self["width"]

    @width.setter
    def width(self, val):
        self["width"] = val

    @property
    def xaxis(self):
        """
        The 'xaxis' property is an instance of XAxis
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.XAxis`
          - A dict of string/value properties that will be passed
            to the XAxis constructor

        Returns
        -------
        plotly.graph_objs.layout.XAxis
        """
        return self["xaxis"]

    @xaxis.setter
    def xaxis(self, val):
        self["xaxis"] = val

    @property
    def yaxis(self):
        """
        The 'yaxis' property is an instance of YAxis
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.YAxis`
          - A dict of string/value properties that will be passed
            to the YAxis constructor

        Returns
        -------
        plotly.graph_objs.layout.YAxis
        """
        return self["yaxis"]

    @yaxis.setter
    def yaxis(self, val):
        self["yaxis"] = val

    @property
    def _prop_descriptions(self):
        return """\
        activeselection
            :class:`plotly.graph_objects.layout.Activeselection`
            instance or dict with compatible properties
        activeshape
            :class:`plotly.graph_objects.layout.Activeshape`
            instance or dict with compatible properties
        annotations
            A tuple of
            :class:`plotly.graph_objects.layout.Annotation`
            instances or dicts with compatible properties
        annotationdefaults
            When used in a template (as
            layout.template.layout.annotationdefaults), sets the
            default property values to use for elements of
            layout.annotations
        autosize
            Determines whether or not a layout width or height that
            has been left undefined by the user is initialized on
            each relayout. Note that, regardless of this attribute,
            an undefined layout width or height is always
            initialized on the first call to plot.
        autotypenumbers
            Using "strict" a numeric string in trace data is not
            converted to a number. Using *convert types* a numeric
            string in trace data may be treated as a number during
            automatic axis `type` detection. This is the default
            value; however it could be overridden for individual
            axes.
        barcornerradius
            Sets the rounding of bar corners. May be an integer
            number of pixels, or a percentage of bar width (as a
            string ending in %).
        bargap
            Sets the gap (in plot fraction) between bars of
            adjacent location coordinates.
        bargroupgap
            Sets the gap (in plot fraction) between bars of the
            same location coordinate.
        barmode
            Determines how bars at the same location coordinate are
            displayed on the graph. With "stack", the bars are
            stacked on top of one another With "relative", the bars
            are stacked on top of one another, with negative values
            below the axis, positive values above With "group", the
            bars are plotted next to one another centered around
            the shared location. With "overlay", the bars are
            plotted over one another, you might need to reduce
            "opacity" to see multiple bars.
        barnorm
            Sets the normalization for bar traces on the graph.
            With "fraction", the value of each bar is divided by
            the sum of all values at that location coordinate.
            "percent" is the same but multiplied by 100 to show
            percentages.
        boxgap
            Sets the gap (in plot fraction) between boxes of
            adjacent location coordinates. Has no effect on traces
            that have "width" set.
        boxgroupgap
            Sets the gap (in plot fraction) between boxes of the
            same location coordinate. Has no effect on traces that
            have "width" set.
        boxmode
            Determines how boxes at the same location coordinate
            are displayed on the graph. If "group", the boxes are
            plotted next to one another centered around the shared
            location. If "overlay", the boxes are plotted over one
            another, you might need to set "opacity" to see them
            multiple boxes. Has no effect on traces that have
            "width" set.
        calendar
            Sets the default calendar system to use for
            interpreting and displaying dates throughout the plot.
        clickmode
            Determines the mode of single click interactions.
            "event" is the default value and emits the
            `plotly_click` event. In addition this mode emits the
            `plotly_selected` event in drag modes "lasso" and
            "select", but with no event data attached (kept for
            compatibility reasons). The "select" flag enables
            selecting single data points via click. This mode also
            supports persistent selections, meaning that pressing
            Shift while clicking, adds to / subtracts from an
            existing selection. "select" with `hovermode`: "x" can
            be confusing, consider explicitly setting `hovermode`:
            "closest" when using this feature. Selection events are
            sent accordingly as long as "event" flag is set as
            well. When the "event" flag is missing, `plotly_click`
            and `plotly_selected` events are not fired.
        coloraxis
            :class:`plotly.graph_objects.layout.Coloraxis` instance
            or dict with compatible properties
        colorscale
            :class:`plotly.graph_objects.layout.Colorscale`
            instance or dict with compatible properties
        colorway
            Sets the default trace colors.
        computed
            Placeholder for exporting automargin-impacting values
            namely `margin.t`, `margin.b`, `margin.l` and
            `margin.r` in "full-json" mode.
        datarevision
            If provided, a changed value tells `Plotly.react` that
            one or more data arrays has changed. This way you can
            modify arrays in-place rather than making a complete
            new copy for an incremental change. If NOT provided,
            `Plotly.react` assumes that data arrays are being
            treated as immutable, thus any data array with a
            different identity from its predecessor contains new
            data.
        dragmode
            Determines the mode of drag interactions. "select" and
            "lasso" apply only to scatter traces with markers or
            text. "orbit" and "turntable" apply only to 3D scenes.
        editrevision
            Controls persistence of user-driven changes in
            `editable: true` configuration, other than trace names
            and axis titles. Defaults to `layout.uirevision`.
        extendfunnelareacolors
            If `true`, the funnelarea slice colors (whether given
            by `funnelareacolorway` or inherited from `colorway`)
            will be extended to three times its original length by
            first repeating every color 20% lighter then each color
            20% darker. This is intended to reduce the likelihood
            of reusing the same color when you have many slices,
            but you can set `false` to disable. Colors provided in
            the trace, using `marker.colors`, are never extended.
        extendiciclecolors
            If `true`, the icicle slice colors (whether given by
            `iciclecolorway` or inherited from `colorway`) will be
            extended to three times its original length by first
            repeating every color 20% lighter then each color 20%
            darker. This is intended to reduce the likelihood of
            reusing the same color when you have many slices, but
            you can set `false` to disable. Colors provided in the
            trace, using `marker.colors`, are never extended.
        extendpiecolors
            If `true`, the pie slice colors (whether given by
            `piecolorway` or inherited from `colorway`) will be
            extended to three times its original length by first
            repeating every color 20% lighter then each color 20%
            darker. This is intended to reduce the likelihood of
            reusing the same color when you have many slices, but
            you can set `false` to disable. Colors provided in the
            trace, using `marker.colors`, are never extended.
        extendsunburstcolors
            If `true`, the sunburst slice colors (whether given by
            `sunburstcolorway` or inherited from `colorway`) will
            be extended to three times its original length by first
            repeating every color 20% lighter then each color 20%
            darker. This is intended to reduce the likelihood of
            reusing the same color when you have many slices, but
            you can set `false` to disable. Colors provided in the
            trace, using `marker.colors`, are never extended.
        extendtreemapcolors
            If `true`, the treemap slice colors (whether given by
            `treemapcolorway` or inherited from `colorway`) will be
            extended to three times its original length by first
            repeating every color 20% lighter then each color 20%
            darker. This is intended to reduce the likelihood of
            reusing the same color when you have many slices, but
            you can set `false` to disable. Colors provided in the
            trace, using `marker.colors`, are never extended.
        font
            Sets the global font. Note that fonts used in traces
            and other layout components inherit from the global
            font.
        funnelareacolorway
            Sets the default funnelarea slice colors. Defaults to
            the main `colorway` used for trace colors. If you
            specify a new list here it can still be extended with
            lighter and darker colors, see
            `extendfunnelareacolors`.
        funnelgap
            Sets the gap (in plot fraction) between bars of
            adjacent location coordinates.
        funnelgroupgap
            Sets the gap (in plot fraction) between bars of the
            same location coordinate.
        funnelmode
            Determines how bars at the same location coordinate are
            displayed on the graph. With "stack", the bars are
            stacked on top of one another With "group", the bars
            are plotted next to one another centered around the
            shared location. With "overlay", the bars are plotted
            over one another, you might need to reduce "opacity" to
            see multiple bars.
        geo
            :class:`plotly.graph_objects.layout.Geo` instance or
            dict with compatible properties
        grid
            :class:`plotly.graph_objects.layout.Grid` instance or
            dict with compatible properties
        height
            Sets the plot's height (in px).
        hiddenlabels
            hiddenlabels is the funnelarea & pie chart analog of
            visible:'legendonly' but it can contain many labels,
            and can simultaneously hide slices from several
            pies/funnelarea charts
        hiddenlabelssrc
            Sets the source reference on Chart Studio Cloud for
            `hiddenlabels`.
        hidesources
            Determines whether or not a text link citing the data
            source is placed at the bottom-right cored of the
            figure. Has only an effect only on graphs that have
            been generated via forked graphs from the Chart Studio
            Cloud (at https://chart-studio.plotly.com or on-
            premise).
        hoverdistance
            Sets the default distance (in pixels) to look for data
            to add hover labels (-1 means no cutoff, 0 means no
            looking for data). This is only a real distance for
            hovering on point-like objects, like scatter points.
            For area-like objects (bars, scatter fills, etc)
            hovering is on inside the area and off outside, but
            these objects will not supersede hover on point-like
            objects in case of conflict.
        hoverlabel
            :class:`plotly.graph_objects.layout.Hoverlabel`
            instance or dict with compatible properties
        hovermode
            Determines the mode of hover interactions. If
            "closest", a single hoverlabel will appear for the
            "closest" point within the `hoverdistance`. If "x" (or
            "y"), multiple hoverlabels will appear for multiple
            points at the "closest" x- (or y-) coordinate within
            the `hoverdistance`, with the caveat that no more than
            one hoverlabel will appear per trace. If *x unified*
            (or *y unified*), a single hoverlabel will appear
            multiple points at the closest x- (or y-) coordinate
            within the `hoverdistance` with the caveat that no more
            than one hoverlabel will appear per trace. In this
            mode, spikelines are enabled by default perpendicular
            to the specified axis. If false, hover interactions are
            disabled.
        hoversubplots
            Determines expansion of hover effects to other subplots
            If "single" just the axis pair of the primary point is
            included without overlaying subplots. If "overlaying"
            all subplots using the main axis and occupying the same
            space are included. If "axis", also include stacked
            subplots using the same axis when `hovermode` is set to
            "x", *x unified*, "y" or *y unified*.
        iciclecolorway
            Sets the default icicle slice colors. Defaults to the
            main `colorway` used for trace colors. If you specify a
            new list here it can still be extended with lighter and
            darker colors, see `extendiciclecolors`.
        images
            A tuple of :class:`plotly.graph_objects.layout.Image`
            instances or dicts with compatible properties
        imagedefaults
            When used in a template (as
            layout.template.layout.imagedefaults), sets the default
            property values to use for elements of layout.images
        legend
            :class:`plotly.graph_objects.layout.Legend` instance or
            dict with compatible properties
        map
            :class:`plotly.graph_objects.layout.Map` instance or
            dict with compatible properties
        mapbox
            :class:`plotly.graph_objects.layout.Mapbox` instance or
            dict with compatible properties
        margin
            :class:`plotly.graph_objects.layout.Margin` instance or
            dict with compatible properties
        meta
            Assigns extra meta information that can be used in
            various `text` attributes. Attributes such as the
            graph, axis and colorbar `title.text`, annotation
            `text` `trace.name` in legend items, `rangeselector`,
            `updatemenus` and `sliders` `label` text all support
            `meta`. One can access `meta` fields using template
            strings: `%{meta[i]}` where `i` is the index of the
            `meta` item in question. `meta` can also be an object
            for example `{key: value}` which can be accessed
            %{meta[key]}.
        metasrc
            Sets the source reference on Chart Studio Cloud for
            `meta`.
        minreducedheight
            Minimum height of the plot with margin.automargin
            applied (in px)
        minreducedwidth
            Minimum width of the plot with margin.automargin
            applied (in px)
        modebar
            :class:`plotly.graph_objects.layout.Modebar` instance
            or dict with compatible properties
        newselection
            :class:`plotly.graph_objects.layout.Newselection`
            instance or dict with compatible properties
        newshape
            :class:`plotly.graph_objects.layout.Newshape` instance
            or dict with compatible properties
        paper_bgcolor
            Sets the background color of the paper where the graph
            is drawn.
        piecolorway
            Sets the default pie slice colors. Defaults to the main
            `colorway` used for trace colors. If you specify a new
            list here it can still be extended with lighter and
            darker colors, see `extendpiecolors`.
        plot_bgcolor
            Sets the background color of the plotting area in-
            between x and y axes.
        polar
            :class:`plotly.graph_objects.layout.Polar` instance or
            dict with compatible properties
        scattergap
            Sets the gap (in plot fraction) between scatter points
            of adjacent location coordinates. Defaults to `bargap`.
        scattermode
            Determines how scatter points at the same location
            coordinate are displayed on the graph. With "group",
            the scatter points are plotted next to one another
            centered around the shared location. With "overlay",
            the scatter points are plotted over one another, you
            might need to reduce "opacity" to see multiple scatter
            points.
        scene
            :class:`plotly.graph_objects.layout.Scene` instance or
            dict with compatible properties
        selectdirection
            When `dragmode` is set to "select", this limits the
            selection of the drag to horizontal, vertical or
            diagonal. "h" only allows horizontal selection, "v"
            only vertical, "d" only diagonal and "any" sets no
            limit.
        selectionrevision
            Controls persistence of user-driven changes in selected
            points from all traces.
        selections
            A tuple of
            :class:`plotly.graph_objects.layout.Selection`
            instances or dicts with compatible properties
        selectiondefaults
            When used in a template (as
            layout.template.layout.selectiondefaults), sets the
            default property values to use for elements of
            layout.selections
        separators
            Sets the decimal and thousand separators. For example,
            *. * puts a '.' before decimals and a space between
            thousands. In English locales, dflt is ".," but other
            locales may alter this default.
        shapes
            A tuple of :class:`plotly.graph_objects.layout.Shape`
            instances or dicts with compatible properties
        shapedefaults
            When used in a template (as
            layout.template.layout.shapedefaults), sets the default
            property values to use for elements of layout.shapes
        showlegend
            Determines whether or not a legend is drawn. Default is
            `true` if there is a trace to show and any of these: a)
            Two or more traces would by default be shown in the
            legend. b) One pie trace is shown in the legend. c) One
            trace is explicitly given with `showlegend: true`.
        sliders
            A tuple of :class:`plotly.graph_objects.layout.Slider`
            instances or dicts with compatible properties
        sliderdefaults
            When used in a template (as
            layout.template.layout.sliderdefaults), sets the
            default property values to use for elements of
            layout.sliders
        smith
            :class:`plotly.graph_objects.layout.Smith` instance or
            dict with compatible properties
        spikedistance
            Sets the default distance (in pixels) to look for data
            to draw spikelines to (-1 means no cutoff, 0 means no
            looking for data). As with hoverdistance, distance does
            not apply to area-like objects. In addition, some
            objects can be hovered on but will not generate
            spikelines, such as scatter fills.
        sunburstcolorway
            Sets the default sunburst slice colors. Defaults to the
            main `colorway` used for trace colors. If you specify a
            new list here it can still be extended with lighter and
            darker colors, see `extendsunburstcolors`.
        template
            Default attributes to be applied to the plot. This
            should be a dict with format: `{'layout':
            layoutTemplate, 'data': {trace_type: [traceTemplate,
            ...], ...}}` where `layoutTemplate` is a dict matching
            the structure of `figure.layout` and `traceTemplate` is
            a dict matching the structure of the trace with type
            `trace_type` (e.g. 'scatter'). Alternatively, this may
            be specified as an instance of
            plotly.graph_objs.layout.Template.  Trace templates are
            applied cyclically to traces of each type. Container
            arrays (eg `annotations`) have special handling: An
            object ending in `defaults` (eg `annotationdefaults`)
            is applied to each array item. But if an item has a
            `templateitemname` key we look in the template array
            for an item with matching `name` and apply that
            instead. If no matching `name` is found we mark the
            item invisible. Any named template item not referenced
            is appended to the end of the array, so this can be
            used to add a watermark annotation or a logo image, for
            example. To omit one of these items on the plot, make
            an item with matching `templateitemname` and `visible:
            false`.
        ternary
            :class:`plotly.graph_objects.layout.Ternary` instance
            or dict with compatible properties
        title
            :class:`plotly.graph_objects.layout.Title` instance or
            dict with compatible properties
        transition
            Sets transition options used during Plotly.react
            updates.
        treemapcolorway
            Sets the default treemap slice colors. Defaults to the
            main `colorway` used for trace colors. If you specify a
            new list here it can still be extended with lighter and
            darker colors, see `extendtreemapcolors`.
        uirevision
            Used to allow user interactions with the plot to
            persist after `Plotly.react` calls that are unaware of
            these interactions. If `uirevision` is omitted, or if
            it is given and it changed from the previous
            `Plotly.react` call, the exact new figure is used. If
            `uirevision` is truthy and did NOT change, any
            attribute that has been affected by user interactions
            and did not receive a different value in the new figure
            will keep the interaction value. `layout.uirevision`
            attribute serves as the default for `uirevision`
            attributes in various sub-containers. For finer control
            you can set these sub-attributes directly. For example,
            if your app separately controls the data on the x and y
            axes you might set `xaxis.uirevision=*time*` and
            `yaxis.uirevision=*cost*`. Then if only the y data is
            changed, you can update `yaxis.uirevision=*quantity*`
            and the y axis range will reset but the x axis range
            will retain any user-driven zoom.
        uniformtext
            :class:`plotly.graph_objects.layout.Uniformtext`
            instance or dict with compatible properties
        updatemenus
            A tuple of
            :class:`plotly.graph_objects.layout.Updatemenu`
            instances or dicts with compatible properties
        updatemenudefaults
            When used in a template (as
            layout.template.layout.updatemenudefaults), sets the
            default property values to use for elements of
            layout.updatemenus
        violingap
            Sets the gap (in plot fraction) between violins of
            adjacent location coordinates. Has no effect on traces
            that have "width" set.
        violingroupgap
            Sets the gap (in plot fraction) between violins of the
            same location coordinate. Has no effect on traces that
            have "width" set.
        violinmode
            Determines how violins at the same location coordinate
            are displayed on the graph. If "group", the violins are
            plotted next to one another centered around the shared
            location. If "overlay", the violins are plotted over
            one another, you might need to set "opacity" to see
            them multiple violins. Has no effect on traces that
            have "width" set.
        waterfallgap
            Sets the gap (in plot fraction) between bars of
            adjacent location coordinates.
        waterfallgroupgap
            Sets the gap (in plot fraction) between bars of the
            same location coordinate.
        waterfallmode
            Determines how bars at the same location coordinate are
            displayed on the graph. With "group", the bars are
            plotted next to one another centered around the shared
            location. With "overlay", the bars are plotted over one
            another, you might need to reduce "opacity" to see
            multiple bars.
        width
            Sets the plot's width (in px).
        xaxis
            :class:`plotly.graph_objects.layout.XAxis` instance or
            dict with compatible properties
        yaxis
            :class:`plotly.graph_objects.layout.YAxis` instance or
            dict with compatible properties
        """

    def __init__(
        self,
        arg=None,
        activeselection=None,
        activeshape=None,
        annotations=None,
        annotationdefaults=None,
        autosize=None,
        autotypenumbers=None,
        barcornerradius=None,
        bargap=None,
        bargroupgap=None,
        barmode=None,
        barnorm=None,
        boxgap=None,
        boxgroupgap=None,
        boxmode=None,
        calendar=None,
        clickmode=None,
        coloraxis=None,
        colorscale=None,
        colorway=None,
        computed=None,
        datarevision=None,
        dragmode=None,
        editrevision=None,
        extendfunnelareacolors=None,
        extendiciclecolors=None,
        extendpiecolors=None,
        extendsunburstcolors=None,
        extendtreemapcolors=None,
        font=None,
        funnelareacolorway=None,
        funnelgap=None,
        funnelgroupgap=None,
        funnelmode=None,
        geo=None,
        grid=None,
        height=None,
        hiddenlabels=None,
        hiddenlabelssrc=None,
        hidesources=None,
        hoverdistance=None,
        hoverlabel=None,
        hovermode=None,
        hoversubplots=None,
        iciclecolorway=None,
        images=None,
        imagedefaults=None,
        legend=None,
        map=None,
        mapbox=None,
        margin=None,
        meta=None,
        metasrc=None,
        minreducedheight=None,
        minreducedwidth=None,
        modebar=None,
        newselection=None,
        newshape=None,
        paper_bgcolor=None,
        piecolorway=None,
        plot_bgcolor=None,
        polar=None,
        scattergap=None,
        scattermode=None,
        scene=None,
        selectdirection=None,
        selectionrevision=None,
        selections=None,
        selectiondefaults=None,
        separators=None,
        shapes=None,
        shapedefaults=None,
        showlegend=None,
        sliders=None,
        sliderdefaults=None,
        smith=None,
        spikedistance=None,
        sunburstcolorway=None,
        template=None,
        ternary=None,
        title=None,
        transition=None,
        treemapcolorway=None,
        uirevision=None,
        uniformtext=None,
        updatemenus=None,
        updatemenudefaults=None,
        violingap=None,
        violingroupgap=None,
        violinmode=None,
        waterfallgap=None,
        waterfallgroupgap=None,
        waterfallmode=None,
        width=None,
        xaxis=None,
        yaxis=None,
        **kwargs,
    ):
        """
        Construct a new Layout object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of :class:`plotly.graph_objs.Layout`
        activeselection
            :class:`plotly.graph_objects.layout.Activeselection`
            instance or dict with compatible properties
        activeshape
            :class:`plotly.graph_objects.layout.Activeshape`
            instance or dict with compatible properties
        annotations
            A tuple of
            :class:`plotly.graph_objects.layout.Annotation`
            instances or dicts with compatible properties
        annotationdefaults
            When used in a template (as
            layout.template.layout.annotationdefaults), sets the
            default property values to use for elements of
            layout.annotations
        autosize
            Determines whether or not a layout width or height that
            has been left undefined by the user is initialized on
            each relayout. Note that, regardless of this attribute,
            an undefined layout width or height is always
            initialized on the first call to plot.
        autotypenumbers
            Using "strict" a numeric string in trace data is not
            converted to a number. Using *convert types* a numeric
            string in trace data may be treated as a number during
            automatic axis `type` detection. This is the default
            value; however it could be overridden for individual
            axes.
        barcornerradius
            Sets the rounding of bar corners. May be an integer
            number of pixels, or a percentage of bar width (as a
            string ending in %).
        bargap
            Sets the gap (in plot fraction) between bars of
            adjacent location coordinates.
        bargroupgap
            Sets the gap (in plot fraction) between bars of the
            same location coordinate.
        barmode
            Determines how bars at the same location coordinate are
            displayed on the graph. With "stack", the bars are
            stacked on top of one another With "relative", the bars
            are stacked on top of one another, with negative values
            below the axis, positive values above With "group", the
            bars are plotted next to one another centered around
            the shared location. With "overlay", the bars are
            plotted over one another, you might need to reduce
            "opacity" to see multiple bars.
        barnorm
            Sets the normalization for bar traces on the graph.
            With "fraction", the value of each bar is divided by
            the sum of all values at that location coordinate.
            "percent" is the same but multiplied by 100 to show
            percentages.
        boxgap
            Sets the gap (in plot fraction) between boxes of
            adjacent location coordinates. Has no effect on traces
            that have "width" set.
        boxgroupgap
            Sets the gap (in plot fraction) between boxes of the
            same location coordinate. Has no effect on traces that
            have "width" set.
        boxmode
            Determines how boxes at the same location coordinate
            are displayed on the graph. If "group", the boxes are
            plotted next to one another centered around the shared
            location. If "overlay", the boxes are plotted over one
            another, you might need to set "opacity" to see them
            multiple boxes. Has no effect on traces that have
            "width" set.
        calendar
            Sets the default calendar system to use for
            interpreting and displaying dates throughout the plot.
        clickmode
            Determines the mode of single click interactions.
            "event" is the default value and emits the
            `plotly_click` event. In addition this mode emits the
            `plotly_selected` event in drag modes "lasso" and
            "select", but with no event data attached (kept for
            compatibility reasons). The "select" flag enables
            selecting single data points via click. This mode also
            supports persistent selections, meaning that pressing
            Shift while clicking, adds to / subtracts from an
            existing selection. "select" with `hovermode`: "x" can
            be confusing, consider explicitly setting `hovermode`:
            "closest" when using this feature. Selection events are
            sent accordingly as long as "event" flag is set as
            well. When the "event" flag is missing, `plotly_click`
            and `plotly_selected` events are not fired.
        coloraxis
            :class:`plotly.graph_objects.layout.Coloraxis` instance
            or dict with compatible properties
        colorscale
            :class:`plotly.graph_objects.layout.Colorscale`
            instance or dict with compatible properties
        colorway
            Sets the default trace colors.
        computed
            Placeholder for exporting automargin-impacting values
            namely `margin.t`, `margin.b`, `margin.l` and
            `margin.r` in "full-json" mode.
        datarevision
            If provided, a changed value tells `Plotly.react` that
            one or more data arrays has changed. This way you can
            modify arrays in-place rather than making a complete
            new copy for an incremental change. If NOT provided,
            `Plotly.react` assumes that data arrays are being
            treated as immutable, thus any data array with a
            different identity from its predecessor contains new
            data.
        dragmode
            Determines the mode of drag interactions. "select" and
            "lasso" apply only to scatter traces with markers or
            text. "orbit" and "turntable" apply only to 3D scenes.
        editrevision
            Controls persistence of user-driven changes in
            `editable: true` configuration, other than trace names
            and axis titles. Defaults to `layout.uirevision`.
        extendfunnelareacolors
            If `true`, the funnelarea slice colors (whether given
            by `funnelareacolorway` or inherited from `colorway`)
            will be extended to three times its original length by
            first repeating every color 20% lighter then each color
            20% darker. This is intended to reduce the likelihood
            of reusing the same color when you have many slices,
            but you can set `false` to disable. Colors provided in
            the trace, using `marker.colors`, are never extended.
        extendiciclecolors
            If `true`, the icicle slice colors (whether given by
            `iciclecolorway` or inherited from `colorway`) will be
            extended to three times its original length by first
            repeating every color 20% lighter then each color 20%
            darker. This is intended to reduce the likelihood of
            reusing the same color when you have many slices, but
            you can set `false` to disable. Colors provided in the
            trace, using `marker.colors`, are never extended.
        extendpiecolors
            If `true`, the pie slice colors (whether given by
            `piecolorway` or inherited from `colorway`) will be
            extended to three times its original length by first
            repeating every color 20% lighter then each color 20%
            darker. This is intended to reduce the likelihood of
            reusing the same color when you have many slices, but
            you can set `false` to disable. Colors provided in the
            trace, using `marker.colors`, are never extended.
        extendsunburstcolors
            If `true`, the sunburst slice colors (whether given by
            `sunburstcolorway` or inherited from `colorway`) will
            be extended to three times its original length by first
            repeating every color 20% lighter then each color 20%
            darker. This is intended to reduce the likelihood of
            reusing the same color when you have many slices, but
            you can set `false` to disable. Colors provided in the
            trace, using `marker.colors`, are never extended.
        extendtreemapcolors
            If `true`, the treemap slice colors (whether given by
            `treemapcolorway` or inherited from `colorway`) will be
            extended to three times its original length by first
            repeating every color 20% lighter then each color 20%
            darker. This is intended to reduce the likelihood of
            reusing the same color when you have many slices, but
            you can set `false` to disable. Colors provided in the
            trace, using `marker.colors`, are never extended.
        font
            Sets the global font. Note that fonts used in traces
            and other layout components inherit from the global
            font.
        funnelareacolorway
            Sets the default funnelarea slice colors. Defaults to
            the main `colorway` used for trace colors. If you
            specify a new list here it can still be extended with
            lighter and darker colors, see
            `extendfunnelareacolors`.
        funnelgap
            Sets the gap (in plot fraction) between bars of
            adjacent location coordinates.
        funnelgroupgap
            Sets the gap (in plot fraction) between bars of the
            same location coordinate.
        funnelmode
            Determines how bars at the same location coordinate are
            displayed on the graph. With "stack", the bars are
            stacked on top of one another With "group", the bars
            are plotted next to one another centered around the
            shared location. With "overlay", the bars are plotted
            over one another, you might need to reduce "opacity" to
            see multiple bars.
        geo
            :class:`plotly.graph_objects.layout.Geo` instance or
            dict with compatible properties
        grid
            :class:`plotly.graph_objects.layout.Grid` instance or
            dict with compatible properties
        height
            Sets the plot's height (in px).
        hiddenlabels
            hiddenlabels is the funnelarea & pie chart analog of
            visible:'legendonly' but it can contain many labels,
            and can simultaneously hide slices from several
            pies/funnelarea charts
        hiddenlabelssrc
            Sets the source reference on Chart Studio Cloud for
            `hiddenlabels`.
        hidesources
            Determines whether or not a text link citing the data
            source is placed at the bottom-right cored of the
            figure. Has only an effect only on graphs that have
            been generated via forked graphs from the Chart Studio
            Cloud (at https://chart-studio.plotly.com or on-
            premise).
        hoverdistance
            Sets the default distance (in pixels) to look for data
            to add hover labels (-1 means no cutoff, 0 means no
            looking for data). This is only a real distance for
            hovering on point-like objects, like scatter points.
            For area-like objects (bars, scatter fills, etc)
            hovering is on inside the area and off outside, but
            these objects will not supersede hover on point-like
            objects in case of conflict.
        hoverlabel
            :class:`plotly.graph_objects.layout.Hoverlabel`
            instance or dict with compatible properties
        hovermode
            Determines the mode of hover interactions. If
            "closest", a single hoverlabel will appear for the
            "closest" point within the `hoverdistance`. If "x" (or
            "y"), multiple hoverlabels will appear for multiple
            points at the "closest" x- (or y-) coordinate within
            the `hoverdistance`, with the caveat that no more than
            one hoverlabel will appear per trace. If *x unified*
            (or *y unified*), a single hoverlabel will appear
            multiple points at the closest x- (or y-) coordinate
            within the `hoverdistance` with the caveat that no more
            than one hoverlabel will appear per trace. In this
            mode, spikelines are enabled by default perpendicular
            to the specified axis. If false, hover interactions are
            disabled.
        hoversubplots
            Determines expansion of hover effects to other subplots
            If "single" just the axis pair of the primary point is
            included without overlaying subplots. If "overlaying"
            all subplots using the main axis and occupying the same
            space are included. If "axis", also include stacked
            subplots using the same axis when `hovermode` is set to
            "x", *x unified*, "y" or *y unified*.
        iciclecolorway
            Sets the default icicle slice colors. Defaults to the
            main `colorway` used for trace colors. If you specify a
            new list here it can still be extended with lighter and
            darker colors, see `extendiciclecolors`.
        images
            A tuple of :class:`plotly.graph_objects.layout.Image`
            instances or dicts with compatible properties
        imagedefaults
            When used in a template (as
            layout.template.layout.imagedefaults), sets the default
            property values to use for elements of layout.images
        legend
            :class:`plotly.graph_objects.layout.Legend` instance or
            dict with compatible properties
        map
            :class:`plotly.graph_objects.layout.Map` instance or
            dict with compatible properties
        mapbox
            :class:`plotly.graph_objects.layout.Mapbox` instance or
            dict with compatible properties
        margin
            :class:`plotly.graph_objects.layout.Margin` instance or
            dict with compatible properties
        meta
            Assigns extra meta information that can be used in
            various `text` attributes. Attributes such as the
            graph, axis and colorbar `title.text`, annotation
            `text` `trace.name` in legend items, `rangeselector`,
            `updatemenus` and `sliders` `label` text all support
            `meta`. One can access `meta` fields using template
            strings: `%{meta[i]}` where `i` is the index of the
            `meta` item in question. `meta` can also be an object
            for example `{key: value}` which can be accessed
            %{meta[key]}.
        metasrc
            Sets the source reference on Chart Studio Cloud for
            `meta`.
        minreducedheight
            Minimum height of the plot with margin.automargin
            applied (in px)
        minreducedwidth
            Minimum width of the plot with margin.automargin
            applied (in px)
        modebar
            :class:`plotly.graph_objects.layout.Modebar` instance
            or dict with compatible properties
        newselection
            :class:`plotly.graph_objects.layout.Newselection`
            instance or dict with compatible properties
        newshape
            :class:`plotly.graph_objects.layout.Newshape` instance
            or dict with compatible properties
        paper_bgcolor
            Sets the background color of the paper where the graph
            is drawn.
        piecolorway
            Sets the default pie slice colors. Defaults to the main
            `colorway` used for trace colors. If you specify a new
            list here it can still be extended with lighter and
            darker colors, see `extendpiecolors`.
        plot_bgcolor
            Sets the background color of the plotting area in-
            between x and y axes.
        polar
            :class:`plotly.graph_objects.layout.Polar` instance or
            dict with compatible properties
        scattergap
            Sets the gap (in plot fraction) between scatter points
            of adjacent location coordinates. Defaults to `bargap`.
        scattermode
            Determines how scatter points at the same location
            coordinate are displayed on the graph. With "group",
            the scatter points are plotted next to one another
            centered around the shared location. With "overlay",
            the scatter points are plotted over one another, you
            might need to reduce "opacity" to see multiple scatter
            points.
        scene
            :class:`plotly.graph_objects.layout.Scene` instance or
            dict with compatible properties
        selectdirection
            When `dragmode` is set to "select", this limits the
            selection of the drag to horizontal, vertical or
            diagonal. "h" only allows horizontal selection, "v"
            only vertical, "d" only diagonal and "any" sets no
            limit.
        selectionrevision
            Controls persistence of user-driven changes in selected
            points from all traces.
        selections
            A tuple of
            :class:`plotly.graph_objects.layout.Selection`
            instances or dicts with compatible properties
        selectiondefaults
            When used in a template (as
            layout.template.layout.selectiondefaults), sets the
            default property values to use for elements of
            layout.selections
        separators
            Sets the decimal and thousand separators. For example,
            *. * puts a '.' before decimals and a space between
            thousands. In English locales, dflt is ".," but other
            locales may alter this default.
        shapes
            A tuple of :class:`plotly.graph_objects.layout.Shape`
            instances or dicts with compatible properties
        shapedefaults
            When used in a template (as
            layout.template.layout.shapedefaults), sets the default
            property values to use for elements of layout.shapes
        showlegend
            Determines whether or not a legend is drawn. Default is
            `true` if there is a trace to show and any of these: a)
            Two or more traces would by default be shown in the
            legend. b) One pie trace is shown in the legend. c) One
            trace is explicitly given with `showlegend: true`.
        sliders
            A tuple of :class:`plotly.graph_objects.layout.Slider`
            instances or dicts with compatible properties
        sliderdefaults
            When used in a template (as
            layout.template.layout.sliderdefaults), sets the
            default property values to use for elements of
            layout.sliders
        smith
            :class:`plotly.graph_objects.layout.Smith` instance or
            dict with compatible properties
        spikedistance
            Sets the default distance (in pixels) to look for data
            to draw spikelines to (-1 means no cutoff, 0 means no
            looking for data). As with hoverdistance, distance does
            not apply to area-like objects. In addition, some
            objects can be hovered on but will not generate
            spikelines, such as scatter fills.
        sunburstcolorway
            Sets the default sunburst slice colors. Defaults to the
            main `colorway` used for trace colors. If you specify a
            new list here it can still be extended with lighter and
            darker colors, see `extendsunburstcolors`.
        template
            Default attributes to be applied to the plot. This
            should be a dict with format: `{'layout':
            layoutTemplate, 'data': {trace_type: [traceTemplate,
            ...], ...}}` where `layoutTemplate` is a dict matching
            the structure of `figure.layout` and `traceTemplate` is
            a dict matching the structure of the trace with type
            `trace_type` (e.g. 'scatter'). Alternatively, this may
            be specified as an instance of
            plotly.graph_objs.layout.Template.  Trace templates are
            applied cyclically to traces of each type. Container
            arrays (eg `annotations`) have special handling: An
            object ending in `defaults` (eg `annotationdefaults`)
            is applied to each array item. But if an item has a
            `templateitemname` key we look in the template array
            for an item with matching `name` and apply that
            instead. If no matching `name` is found we mark the
            item invisible. Any named template item not referenced
            is appended to the end of the array, so this can be
            used to add a watermark annotation or a logo image, for
            example. To omit one of these items on the plot, make
            an item with matching `templateitemname` and `visible:
            false`.
        ternary
            :class:`plotly.graph_objects.layout.Ternary` instance
            or dict with compatible properties
        title
            :class:`plotly.graph_objects.layout.Title` instance or
            dict with compatible properties
        transition
            Sets transition options used during Plotly.react
            updates.
        treemapcolorway
            Sets the default treemap slice colors. Defaults to the
            main `colorway` used for trace colors. If you specify a
            new list here it can still be extended with lighter and
            darker colors, see `extendtreemapcolors`.
        uirevision
            Used to allow user interactions with the plot to
            persist after `Plotly.react` calls that are unaware of
            these interactions. If `uirevision` is omitted, or if
            it is given and it changed from the previous
            `Plotly.react` call, the exact new figure is used. If
            `uirevision` is truthy and did NOT change, any
            attribute that has been affected by user interactions
            and did not receive a different value in the new figure
            will keep the interaction value. `layout.uirevision`
            attribute serves as the default for `uirevision`
            attributes in various sub-containers. For finer control
            you can set these sub-attributes directly. For example,
            if your app separately controls the data on the x and y
            axes you might set `xaxis.uirevision=*time*` and
            `yaxis.uirevision=*cost*`. Then if only the y data is
            changed, you can update `yaxis.uirevision=*quantity*`
            and the y axis range will reset but the x axis range
            will retain any user-driven zoom.
        uniformtext
            :class:`plotly.graph_objects.layout.Uniformtext`
            instance or dict with compatible properties
        updatemenus
            A tuple of
            :class:`plotly.graph_objects.layout.Updatemenu`
            instances or dicts with compatible properties
        updatemenudefaults
            When used in a template (as
            layout.template.layout.updatemenudefaults), sets the
            default property values to use for elements of
            layout.updatemenus
        violingap
            Sets the gap (in plot fraction) between violins of
            adjacent location coordinates. Has no effect on traces
            that have "width" set.
        violingroupgap
            Sets the gap (in plot fraction) between violins of the
            same location coordinate. Has no effect on traces that
            have "width" set.
        violinmode
            Determines how violins at the same location coordinate
            are displayed on the graph. If "group", the violins are
            plotted next to one another centered around the shared
            location. If "overlay", the violins are plotted over
            one another, you might need to set "opacity" to see
            them multiple violins. Has no effect on traces that
            have "width" set.
        waterfallgap
            Sets the gap (in plot fraction) between bars of
            adjacent location coordinates.
        waterfallgroupgap
            Sets the gap (in plot fraction) between bars of the
            same location coordinate.
        waterfallmode
            Determines how bars at the same location coordinate are
            displayed on the graph. With "group", the bars are
            plotted next to one another centered around the shared
            location. With "overlay", the bars are plotted over one
            another, you might need to reduce "opacity" to see
            multiple bars.
        width
            Sets the plot's width (in px).
        xaxis
            :class:`plotly.graph_objects.layout.XAxis` instance or
            dict with compatible properties
        yaxis
            :class:`plotly.graph_objects.layout.YAxis` instance or
            dict with compatible properties

        Returns
        -------
        Layout
        """
        super().__init__("layout")
        if "_parent" in kwargs:
            self._parent = kwargs["_parent"]
            return

        self._valid_props = {
            "activeselection",
            "activeshape",
            "annotationdefaults",
            "annotations",
            "autosize",
            "autotypenumbers",
            "barcornerradius",
            "bargap",
            "bargroupgap",
            "barmode",
            "barnorm",
            "boxgap",
            "boxgroupgap",
            "boxmode",
            "calendar",
            "clickmode",
            "coloraxis",
            "colorscale",
            "colorway",
            "computed",
            "datarevision",
            "dragmode",
            "editrevision",
            "extendfunnelareacolors",
            "extendiciclecolors",
            "extendpiecolors",
            "extendsunburstcolors",
            "extendtreemapcolors",
            "font",
            "funnelareacolorway",
            "funnelgap",
            "funnelgroupgap",
            "funnelmode",
            "geo",
            "grid",
            "height",
            "hiddenlabels",
            "hiddenlabelssrc",
            "hidesources",
            "hoverdistance",
            "hoverlabel",
            "hovermode",
            "hoversubplots",
            "iciclecolorway",
            "imagedefaults",
            "images",
            "legend",
            "map",
            "mapbox",
            "margin",
            "meta",
            "metasrc",
            "minreducedheight",
            "minreducedwidth",
            "modebar",
            "newselection",
            "newshape",
            "paper_bgcolor",
            "piecolorway",
            "plot_bgcolor",
            "polar",
            "scattergap",
            "scattermode",
            "scene",
            "selectdirection",
            "selectiondefaults",
            "selectionrevision",
            "selections",
            "separators",
            "shapedefaults",
            "shapes",
            "showlegend",
            "sliderdefaults",
            "sliders",
            "smith",
            "spikedistance",
            "sunburstcolorway",
            "template",
            "ternary",
            "title",
            "transition",
            "treemapcolorway",
            "uirevision",
            "uniformtext",
            "updatemenudefaults",
            "updatemenus",
            "violingap",
            "violingroupgap",
            "violinmode",
            "waterfallgap",
            "waterfallgroupgap",
            "waterfallmode",
            "width",
            "xaxis",
            "yaxis",
        }

        if arg is None:
            arg = {}
        elif isinstance(arg, self.__class__):
            arg = arg.to_plotly_json()
        elif isinstance(arg, dict):
            arg = _copy.copy(arg)
        else:
            raise ValueError("""\
The first argument to the plotly.graph_objs.Layout
constructor must be a dict or
an instance of :class:`plotly.graph_objs.Layout`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("activeselection", arg, activeselection)
        self._set_property("activeshape", arg, activeshape)
        self._set_property("annotations", arg, annotations)
        self._set_property("annotationdefaults", arg, annotationdefaults)
        self._set_property("autosize", arg, autosize)
        self._set_property("autotypenumbers", arg, autotypenumbers)
        self._set_property("barcornerradius", arg, barcornerradius)
        self._set_property("bargap", arg, bargap)
        self._set_property("bargroupgap", arg, bargroupgap)
        self._set_property("barmode", arg, barmode)
        self._set_property("barnorm", arg, barnorm)
        self._set_property("boxgap", arg, boxgap)
        self._set_property("boxgroupgap", arg, boxgroupgap)
        self._set_property("boxmode", arg, boxmode)
        self._set_property("calendar", arg, calendar)
        self._set_property("clickmode", arg, clickmode)
        self._set_property("coloraxis", arg, coloraxis)
        self._set_property("colorscale", arg, colorscale)
        self._set_property("colorway", arg, colorway)
        self._set_property("computed", arg, computed)
        self._set_property("datarevision", arg, datarevision)
        self._set_property("dragmode", arg, dragmode)
        self._set_property("editrevision", arg, editrevision)
        self._set_property("extendfunnelareacolors", arg, extendfunnelareacolors)
        self._set_property("extendiciclecolors", arg, extendiciclecolors)
        self._set_property("extendpiecolors", arg, extendpiecolors)
        self._set_property("extendsunburstcolors", arg, extendsunburstcolors)
        self._set_property("extendtreemapcolors", arg, extendtreemapcolors)
        self._set_property("font", arg, font)
        self._set_property("funnelareacolorway", arg, funnelareacolorway)
        self._set_property("funnelgap", arg, funnelgap)
        self._set_property("funnelgroupgap", arg, funnelgroupgap)
        self._set_property("funnelmode", arg, funnelmode)
        self._set_property("geo", arg, geo)
        self._set_property("grid", arg, grid)
        self._set_property("height", arg, height)
        self._set_property("hiddenlabels", arg, hiddenlabels)
        self._set_property("hiddenlabelssrc", arg, hiddenlabelssrc)
        self._set_property("hidesources", arg, hidesources)
        self._set_property("hoverdistance", arg, hoverdistance)
        self._set_property("hoverlabel", arg, hoverlabel)
        self._set_property("hovermode", arg, hovermode)
        self._set_property("hoversubplots", arg, hoversubplots)
        self._set_property("iciclecolorway", arg, iciclecolorway)
        self._set_property("images", arg, images)
        self._set_property("imagedefaults", arg, imagedefaults)
        self._set_property("legend", arg, legend)
        self._set_property("map", arg, map)
        self._set_property("mapbox", arg, mapbox)
        self._set_property("margin", arg, margin)
        self._set_property("meta", arg, meta)
        self._set_property("metasrc", arg, metasrc)
        self._set_property("minreducedheight", arg, minreducedheight)
        self._set_property("minreducedwidth", arg, minreducedwidth)
        self._set_property("modebar", arg, modebar)
        self._set_property("newselection", arg, newselection)
        self._set_property("newshape", arg, newshape)
        self._set_property("paper_bgcolor", arg, paper_bgcolor)
        self._set_property("piecolorway", arg, piecolorway)
        self._set_property("plot_bgcolor", arg, plot_bgcolor)
        self._set_property("polar", arg, polar)
        self._set_property("scattergap", arg, scattergap)
        self._set_property("scattermode", arg, scattermode)
        self._set_property("scene", arg, scene)
        self._set_property("selectdirection", arg, selectdirection)
        self._set_property("selectionrevision", arg, selectionrevision)
        self._set_property("selections", arg, selections)
        self._set_property("selectiondefaults", arg, selectiondefaults)
        self._set_property("separators", arg, separators)
        self._set_property("shapes", arg, shapes)
        self._set_property("shapedefaults", arg, shapedefaults)
        self._set_property("showlegend", arg, showlegend)
        self._set_property("sliders", arg, sliders)
        self._set_property("sliderdefaults", arg, sliderdefaults)
        self._set_property("smith", arg, smith)
        self._set_property("spikedistance", arg, spikedistance)
        self._set_property("sunburstcolorway", arg, sunburstcolorway)
        self._set_property("template", arg, template)
        self._set_property("ternary", arg, ternary)
        self._set_property("title", arg, title)
        self._set_property("transition", arg, transition)
        self._set_property("treemapcolorway", arg, treemapcolorway)
        self._set_property("uirevision", arg, uirevision)
        self._set_property("uniformtext", arg, uniformtext)
        self._set_property("updatemenus", arg, updatemenus)
        self._set_property("updatemenudefaults", arg, updatemenudefaults)
        self._set_property("violingap", arg, violingap)
        self._set_property("violingroupgap", arg, violingroupgap)
        self._set_property("violinmode", arg, violinmode)
        self._set_property("waterfallgap", arg, waterfallgap)
        self._set_property("waterfallgroupgap", arg, waterfallgroupgap)
        self._set_property("waterfallmode", arg, waterfallmode)
        self._set_property("width", arg, width)
        self._set_property("xaxis", arg, xaxis)
        self._set_property("yaxis", arg, yaxis)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
