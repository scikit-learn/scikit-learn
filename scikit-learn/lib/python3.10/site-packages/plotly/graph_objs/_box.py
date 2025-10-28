#                   --- THIS FILE IS AUTO-GENERATED ---
# Modifications will be overwitten the next time code generation run.

from plotly.basedatatypes import BaseTraceType as _BaseTraceType
import copy as _copy


class Box(_BaseTraceType):
    _parent_path_str = ""
    _path_str = "box"
    _valid_props = {
        "alignmentgroup",
        "boxmean",
        "boxpoints",
        "customdata",
        "customdatasrc",
        "dx",
        "dy",
        "fillcolor",
        "hoverinfo",
        "hoverinfosrc",
        "hoverlabel",
        "hoveron",
        "hovertemplate",
        "hovertemplatesrc",
        "hovertext",
        "hovertextsrc",
        "ids",
        "idssrc",
        "jitter",
        "legend",
        "legendgroup",
        "legendgrouptitle",
        "legendrank",
        "legendwidth",
        "line",
        "lowerfence",
        "lowerfencesrc",
        "marker",
        "mean",
        "meansrc",
        "median",
        "mediansrc",
        "meta",
        "metasrc",
        "name",
        "notched",
        "notchspan",
        "notchspansrc",
        "notchwidth",
        "offsetgroup",
        "opacity",
        "orientation",
        "pointpos",
        "q1",
        "q1src",
        "q3",
        "q3src",
        "quartilemethod",
        "sd",
        "sdmultiple",
        "sdsrc",
        "selected",
        "selectedpoints",
        "showlegend",
        "showwhiskers",
        "sizemode",
        "stream",
        "text",
        "textsrc",
        "type",
        "uid",
        "uirevision",
        "unselected",
        "upperfence",
        "upperfencesrc",
        "visible",
        "whiskerwidth",
        "width",
        "x",
        "x0",
        "xaxis",
        "xcalendar",
        "xhoverformat",
        "xperiod",
        "xperiod0",
        "xperiodalignment",
        "xsrc",
        "y",
        "y0",
        "yaxis",
        "ycalendar",
        "yhoverformat",
        "yperiod",
        "yperiod0",
        "yperiodalignment",
        "ysrc",
        "zorder",
    }

    @property
    def alignmentgroup(self):
        """
        Set several traces linked to the same position axis or matching
        axes to the same alignmentgroup. This controls whether bars
        compute their positional range dependently or independently.

        The 'alignmentgroup' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["alignmentgroup"]

    @alignmentgroup.setter
    def alignmentgroup(self, val):
        self["alignmentgroup"] = val

    @property
    def boxmean(self):
        """
        If True, the mean of the box(es)' underlying distribution is
        drawn as a dashed line inside the box(es). If "sd" the standard
        deviation is also drawn. Defaults to True when `mean` is set.
        Defaults to "sd" when `sd` is set Otherwise defaults to False.

        The 'boxmean' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                [True, 'sd', False]

        Returns
        -------
        Any
        """
        return self["boxmean"]

    @boxmean.setter
    def boxmean(self, val):
        self["boxmean"] = val

    @property
    def boxpoints(self):
        """
        If "outliers", only the sample points lying outside the
        whiskers are shown If "suspectedoutliers", the outlier points
        are shown and points either less than 4*Q1-3*Q3 or greater than
        4*Q3-3*Q1 are highlighted (see `outliercolor`) If "all", all
        sample points are shown If False, only the box(es) are shown
        with no sample points Defaults to "suspectedoutliers" when
        `marker.outliercolor` or `marker.line.outliercolor` is set.
        Defaults to "all" under the q1/median/q3 signature. Otherwise
        defaults to "outliers".

        The 'boxpoints' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['all', 'outliers', 'suspectedoutliers', False]

        Returns
        -------
        Any
        """
        return self["boxpoints"]

    @boxpoints.setter
    def boxpoints(self, val):
        self["boxpoints"] = val

    @property
    def customdata(self):
        """
        Assigns extra data each datum. This may be useful when
        listening to hover, click and selection events. Note that,
        "scatter" traces also appends customdata items in the markers
        DOM elements

        The 'customdata' property is an array that may be specified as a tuple,
        list, numpy array, or pandas Series

        Returns
        -------
        numpy.ndarray
        """
        return self["customdata"]

    @customdata.setter
    def customdata(self, val):
        self["customdata"] = val

    @property
    def customdatasrc(self):
        """
        Sets the source reference on Chart Studio Cloud for
        `customdata`.

        The 'customdatasrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["customdatasrc"]

    @customdatasrc.setter
    def customdatasrc(self, val):
        self["customdatasrc"] = val

    @property
    def dx(self):
        """
        Sets the x coordinate step for multi-box traces set using
        q1/median/q3.

        The 'dx' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["dx"]

    @dx.setter
    def dx(self, val):
        self["dx"] = val

    @property
    def dy(self):
        """
        Sets the y coordinate step for multi-box traces set using
        q1/median/q3.

        The 'dy' property is a number and may be specified as:
          - An int or float

        Returns
        -------
        int|float
        """
        return self["dy"]

    @dy.setter
    def dy(self, val):
        self["dy"] = val

    @property
    def fillcolor(self):
        """
        Sets the fill color. Defaults to a half-transparent variant of
        the line color, marker color, or marker line color, whichever
        is available.

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
    def hoverinfo(self):
        """
        Determines which trace information appear on hover. If `none`
        or `skip` are set, no information is displayed upon hovering.
        But, if `none` is set, click and hover events are still fired.

        The 'hoverinfo' property is a flaglist and may be specified
        as a string containing:
          - Any combination of ['x', 'y', 'z', 'text', 'name'] joined with '+' characters
            (e.g. 'x+y')
            OR exactly one of ['all', 'none', 'skip'] (e.g. 'skip')
          - A list or array of the above

        Returns
        -------
        Any|numpy.ndarray
        """
        return self["hoverinfo"]

    @hoverinfo.setter
    def hoverinfo(self, val):
        self["hoverinfo"] = val

    @property
    def hoverinfosrc(self):
        """
        Sets the source reference on Chart Studio Cloud for
        `hoverinfo`.

        The 'hoverinfosrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["hoverinfosrc"]

    @hoverinfosrc.setter
    def hoverinfosrc(self, val):
        self["hoverinfosrc"] = val

    @property
    def hoverlabel(self):
        """
        The 'hoverlabel' property is an instance of Hoverlabel
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.box.Hoverlabel`
          - A dict of string/value properties that will be passed
            to the Hoverlabel constructor

        Returns
        -------
        plotly.graph_objs.box.Hoverlabel
        """
        return self["hoverlabel"]

    @hoverlabel.setter
    def hoverlabel(self, val):
        self["hoverlabel"] = val

    @property
    def hoveron(self):
        """
        Do the hover effects highlight individual boxes  or sample
        points or both?

        The 'hoveron' property is a flaglist and may be specified
        as a string containing:
          - Any combination of ['boxes', 'points'] joined with '+' characters
            (e.g. 'boxes+points')

        Returns
        -------
        Any
        """
        return self["hoveron"]

    @hoveron.setter
    def hoveron(self, val):
        self["hoveron"] = val

    @property
    def hovertemplate(self):
        """
        Template string used for rendering the information that appear
        on hover box. Note that this will override `hoverinfo`.
        Variables are inserted using %{variable}, for example "y: %{y}"
        as well as %{xother}, {%_xother}, {%_xother_}, {%xother_}. When
        showing info for several points, "xother" will be added to
        those with different x positions from the first point. An
        underscore before or after "(x|y)other" will add a space on
        that side, only when this field is shown. Numbers are formatted
        using d3-format's syntax %{variable:d3-format}, for example
        "Price: %{y:$.2f}".
        https://github.com/d3/d3-format/tree/v1.4.5#d3-format for
        details on the formatting syntax. Dates are formatted using
        d3-time-format's syntax %{variable|d3-time-format}, for example
        "Day: %{2019-01-01|%A}". https://github.com/d3/d3-time-
        format/tree/v2.2.3#locale_format for details on the date
        formatting syntax. The variables available in `hovertemplate`
        are the ones emitted as event data described at this link
        https://plotly.com/javascript/plotlyjs-events/#event-data.
        Additionally, every attributes that can be specified per-point
        (the ones that are `arrayOk: true`) are available.  Anything
        contained in tag `<extra>` is displayed in the secondary box,
        for example `<extra>%{fullData.name}</extra>`. To hide the
        secondary box completely, use an empty tag `<extra></extra>`.

        The 'hovertemplate' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string
          - A tuple, list, or one-dimensional numpy array of the above

        Returns
        -------
        str|numpy.ndarray
        """
        return self["hovertemplate"]

    @hovertemplate.setter
    def hovertemplate(self, val):
        self["hovertemplate"] = val

    @property
    def hovertemplatesrc(self):
        """
        Sets the source reference on Chart Studio Cloud for
        `hovertemplate`.

        The 'hovertemplatesrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["hovertemplatesrc"]

    @hovertemplatesrc.setter
    def hovertemplatesrc(self, val):
        self["hovertemplatesrc"] = val

    @property
    def hovertext(self):
        """
        Same as `text`.

        The 'hovertext' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string
          - A tuple, list, or one-dimensional numpy array of the above

        Returns
        -------
        str|numpy.ndarray
        """
        return self["hovertext"]

    @hovertext.setter
    def hovertext(self, val):
        self["hovertext"] = val

    @property
    def hovertextsrc(self):
        """
        Sets the source reference on Chart Studio Cloud for
        `hovertext`.

        The 'hovertextsrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["hovertextsrc"]

    @hovertextsrc.setter
    def hovertextsrc(self, val):
        self["hovertextsrc"] = val

    @property
    def ids(self):
        """
        Assigns id labels to each datum. These ids for object constancy
        of data points during animation. Should be an array of strings,
        not numbers or any other type.

        The 'ids' property is an array that may be specified as a tuple,
        list, numpy array, or pandas Series

        Returns
        -------
        numpy.ndarray
        """
        return self["ids"]

    @ids.setter
    def ids(self, val):
        self["ids"] = val

    @property
    def idssrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `ids`.

        The 'idssrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["idssrc"]

    @idssrc.setter
    def idssrc(self, val):
        self["idssrc"] = val

    @property
    def jitter(self):
        """
        Sets the amount of jitter in the sample points drawn. If 0, the
        sample points align along the distribution axis. If 1, the
        sample points are drawn in a random jitter of width equal to
        the width of the box(es).

        The 'jitter' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["jitter"]

    @jitter.setter
    def jitter(self, val):
        self["jitter"] = val

    @property
    def legend(self):
        """
        Sets the reference to a legend to show this trace in.
        References to these legends are "legend", "legend2", "legend3",
        etc. Settings for these legends are set in the layout, under
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
        Sets the legend group for this trace. Traces and shapes part of
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
          - An instance of :class:`plotly.graph_objs.box.Legendgrouptitle`
          - A dict of string/value properties that will be passed
            to the Legendgrouptitle constructor

        Returns
        -------
        plotly.graph_objs.box.Legendgrouptitle
        """
        return self["legendgrouptitle"]

    @legendgrouptitle.setter
    def legendgrouptitle(self, val):
        self["legendgrouptitle"] = val

    @property
    def legendrank(self):
        """
        Sets the legend rank for this trace. Items and groups with
        smaller ranks are presented on top/left side while with
        "reversed" `legend.traceorder` they are on bottom/right side.
        The default legendrank is 1000, so that you can use ranks less
        than 1000 to place certain items before all unranked items, and
        ranks greater than 1000 to go after all unranked items. When
        having unranked or equal rank items shapes would be displayed
        after traces i.e. according to their order in data and layout.

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
        Sets the width (in px or fraction) of the legend for this
        trace.

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
          - An instance of :class:`plotly.graph_objs.box.Line`
          - A dict of string/value properties that will be passed
            to the Line constructor

        Returns
        -------
        plotly.graph_objs.box.Line
        """
        return self["line"]

    @line.setter
    def line(self, val):
        self["line"] = val

    @property
    def lowerfence(self):
        """
        Sets the lower fence values. There should be as many items as
        the number of boxes desired. This attribute has effect only
        under the q1/median/q3 signature. If `lowerfence` is not
        provided but a sample (in `y` or `x`) is set, we compute the
        lower as the last sample point below 1.5 times the IQR.

        The 'lowerfence' property is an array that may be specified as a tuple,
        list, numpy array, or pandas Series

        Returns
        -------
        numpy.ndarray
        """
        return self["lowerfence"]

    @lowerfence.setter
    def lowerfence(self, val):
        self["lowerfence"] = val

    @property
    def lowerfencesrc(self):
        """
        Sets the source reference on Chart Studio Cloud for
        `lowerfence`.

        The 'lowerfencesrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["lowerfencesrc"]

    @lowerfencesrc.setter
    def lowerfencesrc(self, val):
        self["lowerfencesrc"] = val

    @property
    def marker(self):
        """
        The 'marker' property is an instance of Marker
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.box.Marker`
          - A dict of string/value properties that will be passed
            to the Marker constructor

        Returns
        -------
        plotly.graph_objs.box.Marker
        """
        return self["marker"]

    @marker.setter
    def marker(self, val):
        self["marker"] = val

    @property
    def mean(self):
        """
        Sets the mean values. There should be as many items as the
        number of boxes desired. This attribute has effect only under
        the q1/median/q3 signature. If `mean` is not provided but a
        sample (in `y` or `x`) is set, we compute the mean for each box
        using the sample values.

        The 'mean' property is an array that may be specified as a tuple,
        list, numpy array, or pandas Series

        Returns
        -------
        numpy.ndarray
        """
        return self["mean"]

    @mean.setter
    def mean(self, val):
        self["mean"] = val

    @property
    def meansrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `mean`.

        The 'meansrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["meansrc"]

    @meansrc.setter
    def meansrc(self, val):
        self["meansrc"] = val

    @property
    def median(self):
        """
        Sets the median values. There should be as many items as the
        number of boxes desired.

        The 'median' property is an array that may be specified as a tuple,
        list, numpy array, or pandas Series

        Returns
        -------
        numpy.ndarray
        """
        return self["median"]

    @median.setter
    def median(self, val):
        self["median"] = val

    @property
    def mediansrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `median`.

        The 'mediansrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["mediansrc"]

    @mediansrc.setter
    def mediansrc(self, val):
        self["mediansrc"] = val

    @property
    def meta(self):
        """
        Assigns extra meta information associated with this trace that
        can be used in various text attributes. Attributes such as
        trace `name`, graph, axis and colorbar `title.text`, annotation
        `text` `rangeselector`, `updatemenues` and `sliders` `label`
        text all support `meta`. To access the trace `meta` values in
        an attribute in the same trace, simply use `%{meta[i]}` where
        `i` is the index or key of the `meta` item in question. To
        access trace `meta` in layout attributes, use
        `%{data[n[.meta[i]}` where `i` is the index or key of the
        `meta` and `n` is the trace index.

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
    def name(self):
        """
        Sets the trace name. The trace name appears as the legend item
        and on hover. For box traces, the name will also be used for
        the position coordinate, if `x` and `x0` (`y` and `y0` if
        horizontal) are missing and the position axis is categorical

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
    def notched(self):
        """
        Determines whether or not notches are drawn. Notches displays a
        confidence interval around the median. We compute the
        confidence interval as median +/- 1.57 * IQR / sqrt(N), where
        IQR is the interquartile range and N is the sample size. If two
        boxes' notches do not overlap there is 95% confidence their
        medians differ. See
        https://sites.google.com/site/davidsstatistics/home/notched-
        box-plots for more info. Defaults to False unless `notchwidth`
        or `notchspan` is set.

        The 'notched' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["notched"]

    @notched.setter
    def notched(self, val):
        self["notched"] = val

    @property
    def notchspan(self):
        """
        Sets the notch span from the boxes' `median` values. There
        should be as many items as the number of boxes desired. This
        attribute has effect only under the q1/median/q3 signature. If
        `notchspan` is not provided but a sample (in `y` or `x`) is
        set, we compute it as 1.57 * IQR / sqrt(N), where N is the
        sample size.

        The 'notchspan' property is an array that may be specified as a tuple,
        list, numpy array, or pandas Series

        Returns
        -------
        numpy.ndarray
        """
        return self["notchspan"]

    @notchspan.setter
    def notchspan(self, val):
        self["notchspan"] = val

    @property
    def notchspansrc(self):
        """
        Sets the source reference on Chart Studio Cloud for
        `notchspan`.

        The 'notchspansrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["notchspansrc"]

    @notchspansrc.setter
    def notchspansrc(self, val):
        self["notchspansrc"] = val

    @property
    def notchwidth(self):
        """
        Sets the width of the notches relative to the box' width. For
        example, with 0, the notches are as wide as the box(es).

        The 'notchwidth' property is a number and may be specified as:
          - An int or float in the interval [0, 0.5]

        Returns
        -------
        int|float
        """
        return self["notchwidth"]

    @notchwidth.setter
    def notchwidth(self, val):
        self["notchwidth"] = val

    @property
    def offsetgroup(self):
        """
        Set several traces linked to the same position axis or matching
        axes to the same offsetgroup where bars of the same position
        coordinate will line up.

        The 'offsetgroup' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["offsetgroup"]

    @offsetgroup.setter
    def offsetgroup(self, val):
        self["offsetgroup"] = val

    @property
    def opacity(self):
        """
        Sets the opacity of the trace.

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
    def orientation(self):
        """
        Sets the orientation of the box(es). If "v" ("h"), the
        distribution is visualized along the vertical (horizontal).

        The 'orientation' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['v', 'h']

        Returns
        -------
        Any
        """
        return self["orientation"]

    @orientation.setter
    def orientation(self, val):
        self["orientation"] = val

    @property
    def pointpos(self):
        """
        Sets the position of the sample points in relation to the
        box(es). If 0, the sample points are places over the center of
        the box(es). Positive (negative) values correspond to positions
        to the right (left) for vertical boxes and above (below) for
        horizontal boxes

        The 'pointpos' property is a number and may be specified as:
          - An int or float in the interval [-2, 2]

        Returns
        -------
        int|float
        """
        return self["pointpos"]

    @pointpos.setter
    def pointpos(self, val):
        self["pointpos"] = val

    @property
    def q1(self):
        """
        Sets the Quartile 1 values. There should be as many items as
        the number of boxes desired.

        The 'q1' property is an array that may be specified as a tuple,
        list, numpy array, or pandas Series

        Returns
        -------
        numpy.ndarray
        """
        return self["q1"]

    @q1.setter
    def q1(self, val):
        self["q1"] = val

    @property
    def q1src(self):
        """
        Sets the source reference on Chart Studio Cloud for `q1`.

        The 'q1src' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["q1src"]

    @q1src.setter
    def q1src(self, val):
        self["q1src"] = val

    @property
    def q3(self):
        """
        Sets the Quartile 3 values. There should be as many items as
        the number of boxes desired.

        The 'q3' property is an array that may be specified as a tuple,
        list, numpy array, or pandas Series

        Returns
        -------
        numpy.ndarray
        """
        return self["q3"]

    @q3.setter
    def q3(self, val):
        self["q3"] = val

    @property
    def q3src(self):
        """
        Sets the source reference on Chart Studio Cloud for `q3`.

        The 'q3src' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["q3src"]

    @q3src.setter
    def q3src(self, val):
        self["q3src"] = val

    @property
    def quartilemethod(self):
        """
        Sets the method used to compute the sample's Q1 and Q3
        quartiles. The "linear" method uses the 25th percentile for Q1
        and 75th percentile for Q3 as computed using method #10 (listed
        on http://jse.amstat.org/v14n3/langford.html). The "exclusive"
        method uses the median to divide the ordered dataset into two
        halves if the sample is odd, it does not include the median in
        either half - Q1 is then the median of the lower half and Q3
        the median of the upper half. The "inclusive" method also uses
        the median to divide the ordered dataset into two halves but if
        the sample is odd, it includes the median in both halves - Q1
        is then the median of the lower half and Q3 the median of the
        upper half.

        The 'quartilemethod' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['linear', 'exclusive', 'inclusive']

        Returns
        -------
        Any
        """
        return self["quartilemethod"]

    @quartilemethod.setter
    def quartilemethod(self, val):
        self["quartilemethod"] = val

    @property
    def sd(self):
        """
        Sets the standard deviation values. There should be as many
        items as the number of boxes desired. This attribute has effect
        only under the q1/median/q3 signature. If `sd` is not provided
        but a sample (in `y` or `x`) is set, we compute the standard
        deviation for each box using the sample values.

        The 'sd' property is an array that may be specified as a tuple,
        list, numpy array, or pandas Series

        Returns
        -------
        numpy.ndarray
        """
        return self["sd"]

    @sd.setter
    def sd(self, val):
        self["sd"] = val

    @property
    def sdmultiple(self):
        """
        Scales the box size when sizemode=sd Allowing boxes to be drawn
        across any stddev range For example 1-stddev, 3-stddev,
        5-stddev

        The 'sdmultiple' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["sdmultiple"]

    @sdmultiple.setter
    def sdmultiple(self, val):
        self["sdmultiple"] = val

    @property
    def sdsrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `sd`.

        The 'sdsrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["sdsrc"]

    @sdsrc.setter
    def sdsrc(self, val):
        self["sdsrc"] = val

    @property
    def selected(self):
        """
        The 'selected' property is an instance of Selected
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.box.Selected`
          - A dict of string/value properties that will be passed
            to the Selected constructor

        Returns
        -------
        plotly.graph_objs.box.Selected
        """
        return self["selected"]

    @selected.setter
    def selected(self, val):
        self["selected"] = val

    @property
    def selectedpoints(self):
        """
        Array containing integer indices of selected points. Has an
        effect only for traces that support selections. Note that an
        empty array means an empty selection where the `unselected` are
        turned on for all points, whereas, any other non-array values
        means no selection all where the `selected` and `unselected`
        styles have no effect.

        The 'selectedpoints' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["selectedpoints"]

    @selectedpoints.setter
    def selectedpoints(self, val):
        self["selectedpoints"] = val

    @property
    def showlegend(self):
        """
        Determines whether or not an item corresponding to this trace
        is shown in the legend.

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
    def showwhiskers(self):
        """
        Determines whether or not whiskers are visible. Defaults to
        true for `sizemode` "quartiles", false for "sd".

        The 'showwhiskers' property must be specified as a bool
        (either True, or False)

        Returns
        -------
        bool
        """
        return self["showwhiskers"]

    @showwhiskers.setter
    def showwhiskers(self, val):
        self["showwhiskers"] = val

    @property
    def sizemode(self):
        """
        Sets the upper and lower bound for the boxes quartiles means
        box is drawn between Q1 and Q3 SD means the box is drawn
        between Mean +- Standard Deviation Argument sdmultiple (default
        1) to scale the box size So it could be drawn 1-stddev,
        3-stddev etc

        The 'sizemode' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['quartiles', 'sd']

        Returns
        -------
        Any
        """
        return self["sizemode"]

    @sizemode.setter
    def sizemode(self, val):
        self["sizemode"] = val

    @property
    def stream(self):
        """
        The 'stream' property is an instance of Stream
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.box.Stream`
          - A dict of string/value properties that will be passed
            to the Stream constructor

        Returns
        -------
        plotly.graph_objs.box.Stream
        """
        return self["stream"]

    @stream.setter
    def stream(self, val):
        self["stream"] = val

    @property
    def text(self):
        """
        Sets the text elements associated with each sample value. If a
        single string, the same string appears over all the data
        points. If an array of string, the items are mapped in order to
        the this trace's (x,y) coordinates. To be seen, trace
        `hoverinfo` must contain a "text" flag.

        The 'text' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string
          - A tuple, list, or one-dimensional numpy array of the above

        Returns
        -------
        str|numpy.ndarray
        """
        return self["text"]

    @text.setter
    def text(self, val):
        self["text"] = val

    @property
    def textsrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `text`.

        The 'textsrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["textsrc"]

    @textsrc.setter
    def textsrc(self, val):
        self["textsrc"] = val

    @property
    def uid(self):
        """
        Assign an id to this trace, Use this to provide object
        constancy between traces during animations and transitions.

        The 'uid' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["uid"]

    @uid.setter
    def uid(self, val):
        self["uid"] = val

    @property
    def uirevision(self):
        """
        Controls persistence of some user-driven changes to the trace:
        `constraintrange` in `parcoords` traces, as well as some
        `editable: true` modifications such as `name` and
        `colorbar.title`. Defaults to `layout.uirevision`. Note that
        other user-driven trace attribute changes are controlled by
        `layout` attributes: `trace.visible` is controlled by
        `layout.legend.uirevision`, `selectedpoints` is controlled by
        `layout.selectionrevision`, and `colorbar.(x|y)` (accessible
        with `config: {editable: true}`) is controlled by
        `layout.editrevision`. Trace changes are tracked by `uid`,
        which only falls back on trace index if no `uid` is provided.
        So if your app can add/remove traces before the end of the
        `data` array, such that the same trace has a different index,
        you can still preserve user-driven changes if you give each
        trace a `uid` that stays with it as it moves.

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
    def unselected(self):
        """
        The 'unselected' property is an instance of Unselected
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.box.Unselected`
          - A dict of string/value properties that will be passed
            to the Unselected constructor

        Returns
        -------
        plotly.graph_objs.box.Unselected
        """
        return self["unselected"]

    @unselected.setter
    def unselected(self, val):
        self["unselected"] = val

    @property
    def upperfence(self):
        """
        Sets the upper fence values. There should be as many items as
        the number of boxes desired. This attribute has effect only
        under the q1/median/q3 signature. If `upperfence` is not
        provided but a sample (in `y` or `x`) is set, we compute the
        upper as the last sample point above 1.5 times the IQR.

        The 'upperfence' property is an array that may be specified as a tuple,
        list, numpy array, or pandas Series

        Returns
        -------
        numpy.ndarray
        """
        return self["upperfence"]

    @upperfence.setter
    def upperfence(self, val):
        self["upperfence"] = val

    @property
    def upperfencesrc(self):
        """
        Sets the source reference on Chart Studio Cloud for
        `upperfence`.

        The 'upperfencesrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["upperfencesrc"]

    @upperfencesrc.setter
    def upperfencesrc(self, val):
        self["upperfencesrc"] = val

    @property
    def visible(self):
        """
        Determines whether or not this trace is visible. If
        "legendonly", the trace is not drawn, but can appear as a
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
    def whiskerwidth(self):
        """
        Sets the width of the whiskers relative to the box' width. For
        example, with 1, the whiskers are as wide as the box(es).

        The 'whiskerwidth' property is a number and may be specified as:
          - An int or float in the interval [0, 1]

        Returns
        -------
        int|float
        """
        return self["whiskerwidth"]

    @whiskerwidth.setter
    def whiskerwidth(self, val):
        self["whiskerwidth"] = val

    @property
    def width(self):
        """
        Sets the width of the box in data coordinate If 0 (default
        value) the width is automatically selected based on the
        positions of other box traces in the same subplot.

        The 'width' property is a number and may be specified as:
          - An int or float in the interval [0, inf]

        Returns
        -------
        int|float
        """
        return self["width"]

    @width.setter
    def width(self, val):
        self["width"] = val

    @property
    def x(self):
        """
        Sets the x sample data or coordinates. See overview for more
        info.

        The 'x' property is an array that may be specified as a tuple,
        list, numpy array, or pandas Series

        Returns
        -------
        numpy.ndarray
        """
        return self["x"]

    @x.setter
    def x(self, val):
        self["x"] = val

    @property
    def x0(self):
        """
        Sets the x coordinate for single-box traces or the starting
        coordinate for multi-box traces set using q1/median/q3. See
        overview for more info.

        The 'x0' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["x0"]

    @x0.setter
    def x0(self, val):
        self["x0"] = val

    @property
    def xaxis(self):
        """
        Sets a reference between this trace's x coordinates and a 2D
        cartesian x axis. If "x" (the default value), the x coordinates
        refer to `layout.xaxis`. If "x2", the x coordinates refer to
        `layout.xaxis2`, and so on.

        The 'xaxis' property is an identifier of a particular
        subplot, of type 'x', that may be specified as the string 'x'
        optionally followed by an integer >= 1
        (e.g. 'x', 'x1', 'x2', 'x3', etc.)

        Returns
        -------
        str
        """
        return self["xaxis"]

    @xaxis.setter
    def xaxis(self, val):
        self["xaxis"] = val

    @property
    def xcalendar(self):
        """
        Sets the calendar system to use with `x` date data.

        The 'xcalendar' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['chinese', 'coptic', 'discworld', 'ethiopian',
                'gregorian', 'hebrew', 'islamic', 'jalali', 'julian',
                'mayan', 'nanakshahi', 'nepali', 'persian', 'taiwan',
                'thai', 'ummalqura']

        Returns
        -------
        Any
        """
        return self["xcalendar"]

    @xcalendar.setter
    def xcalendar(self, val):
        self["xcalendar"] = val

    @property
    def xhoverformat(self):
        """
        Sets the hover text formatting rulefor `x`  using d3 formatting
        mini-languages which are very similar to those in Python. For
        numbers, see:
        https://github.com/d3/d3-format/tree/v1.4.5#d3-format. And for
        dates see: https://github.com/d3/d3-time-
        format/tree/v2.2.3#locale_format. We add two items to d3's date
        formatter: "%h" for half of the year as a decimal number as
        well as "%{n}f" for fractional seconds with n digits. For
        example, *2016-10-13 09:15:23.456* with tickformat
        "%H~%M~%S.%2f" would display *09~15~23.46*By default the values
        are formatted using `xaxis.hoverformat`.

        The 'xhoverformat' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["xhoverformat"]

    @xhoverformat.setter
    def xhoverformat(self, val):
        self["xhoverformat"] = val

    @property
    def xperiod(self):
        """
        Only relevant when the axis `type` is "date". Sets the period
        positioning in milliseconds or "M<n>" on the x axis. Special
        values in the form of "M<n>" could be used to declare the
        number of months. In this case `n` must be a positive integer.

        The 'xperiod' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["xperiod"]

    @xperiod.setter
    def xperiod(self, val):
        self["xperiod"] = val

    @property
    def xperiod0(self):
        """
        Only relevant when the axis `type` is "date". Sets the base for
        period positioning in milliseconds or date string on the x0
        axis. When `x0period` is round number of weeks, the `x0period0`
        by default would be on a Sunday i.e. 2000-01-02, otherwise it
        would be at 2000-01-01.

        The 'xperiod0' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["xperiod0"]

    @xperiod0.setter
    def xperiod0(self, val):
        self["xperiod0"] = val

    @property
    def xperiodalignment(self):
        """
        Only relevant when the axis `type` is "date". Sets the
        alignment of data points on the x axis.

        The 'xperiodalignment' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['start', 'middle', 'end']

        Returns
        -------
        Any
        """
        return self["xperiodalignment"]

    @xperiodalignment.setter
    def xperiodalignment(self, val):
        self["xperiodalignment"] = val

    @property
    def xsrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `x`.

        The 'xsrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["xsrc"]

    @xsrc.setter
    def xsrc(self, val):
        self["xsrc"] = val

    @property
    def y(self):
        """
        Sets the y sample data or coordinates. See overview for more
        info.

        The 'y' property is an array that may be specified as a tuple,
        list, numpy array, or pandas Series

        Returns
        -------
        numpy.ndarray
        """
        return self["y"]

    @y.setter
    def y(self, val):
        self["y"] = val

    @property
    def y0(self):
        """
        Sets the y coordinate for single-box traces or the starting
        coordinate for multi-box traces set using q1/median/q3. See
        overview for more info.

        The 'y0' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["y0"]

    @y0.setter
    def y0(self, val):
        self["y0"] = val

    @property
    def yaxis(self):
        """
        Sets a reference between this trace's y coordinates and a 2D
        cartesian y axis. If "y" (the default value), the y coordinates
        refer to `layout.yaxis`. If "y2", the y coordinates refer to
        `layout.yaxis2`, and so on.

        The 'yaxis' property is an identifier of a particular
        subplot, of type 'y', that may be specified as the string 'y'
        optionally followed by an integer >= 1
        (e.g. 'y', 'y1', 'y2', 'y3', etc.)

        Returns
        -------
        str
        """
        return self["yaxis"]

    @yaxis.setter
    def yaxis(self, val):
        self["yaxis"] = val

    @property
    def ycalendar(self):
        """
        Sets the calendar system to use with `y` date data.

        The 'ycalendar' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['chinese', 'coptic', 'discworld', 'ethiopian',
                'gregorian', 'hebrew', 'islamic', 'jalali', 'julian',
                'mayan', 'nanakshahi', 'nepali', 'persian', 'taiwan',
                'thai', 'ummalqura']

        Returns
        -------
        Any
        """
        return self["ycalendar"]

    @ycalendar.setter
    def ycalendar(self, val):
        self["ycalendar"] = val

    @property
    def yhoverformat(self):
        """
        Sets the hover text formatting rulefor `y`  using d3 formatting
        mini-languages which are very similar to those in Python. For
        numbers, see:
        https://github.com/d3/d3-format/tree/v1.4.5#d3-format. And for
        dates see: https://github.com/d3/d3-time-
        format/tree/v2.2.3#locale_format. We add two items to d3's date
        formatter: "%h" for half of the year as a decimal number as
        well as "%{n}f" for fractional seconds with n digits. For
        example, *2016-10-13 09:15:23.456* with tickformat
        "%H~%M~%S.%2f" would display *09~15~23.46*By default the values
        are formatted using `yaxis.hoverformat`.

        The 'yhoverformat' property is a string and must be specified as:
          - A string
          - A number that will be converted to a string

        Returns
        -------
        str
        """
        return self["yhoverformat"]

    @yhoverformat.setter
    def yhoverformat(self, val):
        self["yhoverformat"] = val

    @property
    def yperiod(self):
        """
        Only relevant when the axis `type` is "date". Sets the period
        positioning in milliseconds or "M<n>" on the y axis. Special
        values in the form of "M<n>" could be used to declare the
        number of months. In this case `n` must be a positive integer.

        The 'yperiod' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["yperiod"]

    @yperiod.setter
    def yperiod(self, val):
        self["yperiod"] = val

    @property
    def yperiod0(self):
        """
        Only relevant when the axis `type` is "date". Sets the base for
        period positioning in milliseconds or date string on the y0
        axis. When `y0period` is round number of weeks, the `y0period0`
        by default would be on a Sunday i.e. 2000-01-02, otherwise it
        would be at 2000-01-01.

        The 'yperiod0' property accepts values of any type

        Returns
        -------
        Any
        """
        return self["yperiod0"]

    @yperiod0.setter
    def yperiod0(self, val):
        self["yperiod0"] = val

    @property
    def yperiodalignment(self):
        """
        Only relevant when the axis `type` is "date". Sets the
        alignment of data points on the y axis.

        The 'yperiodalignment' property is an enumeration that may be specified as:
          - One of the following enumeration values:
                ['start', 'middle', 'end']

        Returns
        -------
        Any
        """
        return self["yperiodalignment"]

    @yperiodalignment.setter
    def yperiodalignment(self, val):
        self["yperiodalignment"] = val

    @property
    def ysrc(self):
        """
        Sets the source reference on Chart Studio Cloud for `y`.

        The 'ysrc' property must be specified as a string or
        as a plotly.grid_objs.Column object

        Returns
        -------
        str
        """
        return self["ysrc"]

    @ysrc.setter
    def ysrc(self, val):
        self["ysrc"] = val

    @property
    def zorder(self):
        """
        Sets the layer on which this trace is displayed, relative to
        other SVG traces on the same subplot. SVG traces with higher
        `zorder` appear in front of those with lower `zorder`.

        The 'zorder' property is a integer and may be specified as:
          - An int (or float that will be cast to an int)

        Returns
        -------
        int
        """
        return self["zorder"]

    @zorder.setter
    def zorder(self, val):
        self["zorder"] = val

    @property
    def type(self):
        return self._props["type"]

    @property
    def _prop_descriptions(self):
        return """\
        alignmentgroup
            Set several traces linked to the same position axis or
            matching axes to the same alignmentgroup. This controls
            whether bars compute their positional range dependently
            or independently.
        boxmean
            If True, the mean of the box(es)' underlying
            distribution is drawn as a dashed line inside the
            box(es). If "sd" the standard deviation is also drawn.
            Defaults to True when `mean` is set. Defaults to "sd"
            when `sd` is set Otherwise defaults to False.
        boxpoints
            If "outliers", only the sample points lying outside the
            whiskers are shown If "suspectedoutliers", the outlier
            points are shown and points either less than 4*Q1-3*Q3
            or greater than 4*Q3-3*Q1 are highlighted (see
            `outliercolor`) If "all", all sample points are shown
            If False, only the box(es) are shown with no sample
            points Defaults to "suspectedoutliers" when
            `marker.outliercolor` or `marker.line.outliercolor` is
            set. Defaults to "all" under the q1/median/q3
            signature. Otherwise defaults to "outliers".
        customdata
            Assigns extra data each datum. This may be useful when
            listening to hover, click and selection events. Note
            that, "scatter" traces also appends customdata items in
            the markers DOM elements
        customdatasrc
            Sets the source reference on Chart Studio Cloud for
            `customdata`.
        dx
            Sets the x coordinate step for multi-box traces set
            using q1/median/q3.
        dy
            Sets the y coordinate step for multi-box traces set
            using q1/median/q3.
        fillcolor
            Sets the fill color. Defaults to a half-transparent
            variant of the line color, marker color, or marker line
            color, whichever is available.
        hoverinfo
            Determines which trace information appear on hover. If
            `none` or `skip` are set, no information is displayed
            upon hovering. But, if `none` is set, click and hover
            events are still fired.
        hoverinfosrc
            Sets the source reference on Chart Studio Cloud for
            `hoverinfo`.
        hoverlabel
            :class:`plotly.graph_objects.box.Hoverlabel` instance
            or dict with compatible properties
        hoveron
            Do the hover effects highlight individual boxes  or
            sample points or both?
        hovertemplate
            Template string used for rendering the information that
            appear on hover box. Note that this will override
            `hoverinfo`. Variables are inserted using %{variable},
            for example "y: %{y}" as well as %{xother}, {%_xother},
            {%_xother_}, {%xother_}. When showing info for several
            points, "xother" will be added to those with different
            x positions from the first point. An underscore before
            or after "(x|y)other" will add a space on that side,
            only when this field is shown. Numbers are formatted
            using d3-format's syntax %{variable:d3-format}, for
            example "Price: %{y:$.2f}".
            https://github.com/d3/d3-format/tree/v1.4.5#d3-format
            for details on the formatting syntax. Dates are
            formatted using d3-time-format's syntax
            %{variable|d3-time-format}, for example "Day:
            %{2019-01-01|%A}". https://github.com/d3/d3-time-
            format/tree/v2.2.3#locale_format for details on the
            date formatting syntax. The variables available in
            `hovertemplate` are the ones emitted as event data
            described at this link
            https://plotly.com/javascript/plotlyjs-events/#event-
            data. Additionally, every attributes that can be
            specified per-point (the ones that are `arrayOk: true`)
            are available.  Anything contained in tag `<extra>` is
            displayed in the secondary box, for example
            `<extra>%{fullData.name}</extra>`. To hide the
            secondary box completely, use an empty tag
            `<extra></extra>`.
        hovertemplatesrc
            Sets the source reference on Chart Studio Cloud for
            `hovertemplate`.
        hovertext
            Same as `text`.
        hovertextsrc
            Sets the source reference on Chart Studio Cloud for
            `hovertext`.
        ids
            Assigns id labels to each datum. These ids for object
            constancy of data points during animation. Should be an
            array of strings, not numbers or any other type.
        idssrc
            Sets the source reference on Chart Studio Cloud for
            `ids`.
        jitter
            Sets the amount of jitter in the sample points drawn.
            If 0, the sample points align along the distribution
            axis. If 1, the sample points are drawn in a random
            jitter of width equal to the width of the box(es).
        legend
            Sets the reference to a legend to show this trace in.
            References to these legends are "legend", "legend2",
            "legend3", etc. Settings for these legends are set in
            the layout, under `layout.legend`, `layout.legend2`,
            etc.
        legendgroup
            Sets the legend group for this trace. Traces and shapes
            part of the same legend group hide/show at the same
            time when toggling legend items.
        legendgrouptitle
            :class:`plotly.graph_objects.box.Legendgrouptitle`
            instance or dict with compatible properties
        legendrank
            Sets the legend rank for this trace. Items and groups
            with smaller ranks are presented on top/left side while
            with "reversed" `legend.traceorder` they are on
            bottom/right side. The default legendrank is 1000, so
            that you can use ranks less than 1000 to place certain
            items before all unranked items, and ranks greater than
            1000 to go after all unranked items. When having
            unranked or equal rank items shapes would be displayed
            after traces i.e. according to their order in data and
            layout.
        legendwidth
            Sets the width (in px or fraction) of the legend for
            this trace.
        line
            :class:`plotly.graph_objects.box.Line` instance or dict
            with compatible properties
        lowerfence
            Sets the lower fence values. There should be as many
            items as the number of boxes desired. This attribute
            has effect only under the q1/median/q3 signature. If
            `lowerfence` is not provided but a sample (in `y` or
            `x`) is set, we compute the lower as the last sample
            point below 1.5 times the IQR.
        lowerfencesrc
            Sets the source reference on Chart Studio Cloud for
            `lowerfence`.
        marker
            :class:`plotly.graph_objects.box.Marker` instance or
            dict with compatible properties
        mean
            Sets the mean values. There should be as many items as
            the number of boxes desired. This attribute has effect
            only under the q1/median/q3 signature. If `mean` is not
            provided but a sample (in `y` or `x`) is set, we
            compute the mean for each box using the sample values.
        meansrc
            Sets the source reference on Chart Studio Cloud for
            `mean`.
        median
            Sets the median values. There should be as many items
            as the number of boxes desired.
        mediansrc
            Sets the source reference on Chart Studio Cloud for
            `median`.
        meta
            Assigns extra meta information associated with this
            trace that can be used in various text attributes.
            Attributes such as trace `name`, graph, axis and
            colorbar `title.text`, annotation `text`
            `rangeselector`, `updatemenues` and `sliders` `label`
            text all support `meta`. To access the trace `meta`
            values in an attribute in the same trace, simply use
            `%{meta[i]}` where `i` is the index or key of the
            `meta` item in question. To access trace `meta` in
            layout attributes, use `%{data[n[.meta[i]}` where `i`
            is the index or key of the `meta` and `n` is the trace
            index.
        metasrc
            Sets the source reference on Chart Studio Cloud for
            `meta`.
        name
            Sets the trace name. The trace name appears as the
            legend item and on hover. For box traces, the name will
            also be used for the position coordinate, if `x` and
            `x0` (`y` and `y0` if horizontal) are missing and the
            position axis is categorical
        notched
            Determines whether or not notches are drawn. Notches
            displays a confidence interval around the median. We
            compute the confidence interval as median +/- 1.57 *
            IQR / sqrt(N), where IQR is the interquartile range and
            N is the sample size. If two boxes' notches do not
            overlap there is 95% confidence their medians differ.
            See https://sites.google.com/site/davidsstatistics/home
            /notched-box-plots for more info. Defaults to False
            unless `notchwidth` or `notchspan` is set.
        notchspan
            Sets the notch span from the boxes' `median` values.
            There should be as many items as the number of boxes
            desired. This attribute has effect only under the
            q1/median/q3 signature. If `notchspan` is not provided
            but a sample (in `y` or `x`) is set, we compute it as
            1.57 * IQR / sqrt(N), where N is the sample size.
        notchspansrc
            Sets the source reference on Chart Studio Cloud for
            `notchspan`.
        notchwidth
            Sets the width of the notches relative to the box'
            width. For example, with 0, the notches are as wide as
            the box(es).
        offsetgroup
            Set several traces linked to the same position axis or
            matching axes to the same offsetgroup where bars of the
            same position coordinate will line up.
        opacity
            Sets the opacity of the trace.
        orientation
            Sets the orientation of the box(es). If "v" ("h"), the
            distribution is visualized along the vertical
            (horizontal).
        pointpos
            Sets the position of the sample points in relation to
            the box(es). If 0, the sample points are places over
            the center of the box(es). Positive (negative) values
            correspond to positions to the right (left) for
            vertical boxes and above (below) for horizontal boxes
        q1
            Sets the Quartile 1 values. There should be as many
            items as the number of boxes desired.
        q1src
            Sets the source reference on Chart Studio Cloud for
            `q1`.
        q3
            Sets the Quartile 3 values. There should be as many
            items as the number of boxes desired.
        q3src
            Sets the source reference on Chart Studio Cloud for
            `q3`.
        quartilemethod
            Sets the method used to compute the sample's Q1 and Q3
            quartiles. The "linear" method uses the 25th percentile
            for Q1 and 75th percentile for Q3 as computed using
            method #10 (listed on
            http://jse.amstat.org/v14n3/langford.html). The
            "exclusive" method uses the median to divide the
            ordered dataset into two halves if the sample is odd,
            it does not include the median in either half - Q1 is
            then the median of the lower half and Q3 the median of
            the upper half. The "inclusive" method also uses the
            median to divide the ordered dataset into two halves
            but if the sample is odd, it includes the median in
            both halves - Q1 is then the median of the lower half
            and Q3 the median of the upper half.
        sd
            Sets the standard deviation values. There should be as
            many items as the number of boxes desired. This
            attribute has effect only under the q1/median/q3
            signature. If `sd` is not provided but a sample (in `y`
            or `x`) is set, we compute the standard deviation for
            each box using the sample values.
        sdmultiple
            Scales the box size when sizemode=sd Allowing boxes to
            be drawn across any stddev range For example 1-stddev,
            3-stddev, 5-stddev
        sdsrc
            Sets the source reference on Chart Studio Cloud for
            `sd`.
        selected
            :class:`plotly.graph_objects.box.Selected` instance or
            dict with compatible properties
        selectedpoints
            Array containing integer indices of selected points.
            Has an effect only for traces that support selections.
            Note that an empty array means an empty selection where
            the `unselected` are turned on for all points, whereas,
            any other non-array values means no selection all where
            the `selected` and `unselected` styles have no effect.
        showlegend
            Determines whether or not an item corresponding to this
            trace is shown in the legend.
        showwhiskers
            Determines whether or not whiskers are visible.
            Defaults to true for `sizemode` "quartiles", false for
            "sd".
        sizemode
            Sets the upper and lower bound for the boxes quartiles
            means box is drawn between Q1 and Q3 SD means the box
            is drawn between Mean +- Standard Deviation Argument
            sdmultiple (default 1) to scale the box size So it
            could be drawn 1-stddev, 3-stddev etc
        stream
            :class:`plotly.graph_objects.box.Stream` instance or
            dict with compatible properties
        text
            Sets the text elements associated with each sample
            value. If a single string, the same string appears over
            all the data points. If an array of string, the items
            are mapped in order to the this trace's (x,y)
            coordinates. To be seen, trace `hoverinfo` must contain
            a "text" flag.
        textsrc
            Sets the source reference on Chart Studio Cloud for
            `text`.
        uid
            Assign an id to this trace, Use this to provide object
            constancy between traces during animations and
            transitions.
        uirevision
            Controls persistence of some user-driven changes to the
            trace: `constraintrange` in `parcoords` traces, as well
            as some `editable: true` modifications such as `name`
            and `colorbar.title`. Defaults to `layout.uirevision`.
            Note that other user-driven trace attribute changes are
            controlled by `layout` attributes: `trace.visible` is
            controlled by `layout.legend.uirevision`,
            `selectedpoints` is controlled by
            `layout.selectionrevision`, and `colorbar.(x|y)`
            (accessible with `config: {editable: true}`) is
            controlled by `layout.editrevision`. Trace changes are
            tracked by `uid`, which only falls back on trace index
            if no `uid` is provided. So if your app can add/remove
            traces before the end of the `data` array, such that
            the same trace has a different index, you can still
            preserve user-driven changes if you give each trace a
            `uid` that stays with it as it moves.
        unselected
            :class:`plotly.graph_objects.box.Unselected` instance
            or dict with compatible properties
        upperfence
            Sets the upper fence values. There should be as many
            items as the number of boxes desired. This attribute
            has effect only under the q1/median/q3 signature. If
            `upperfence` is not provided but a sample (in `y` or
            `x`) is set, we compute the upper as the last sample
            point above 1.5 times the IQR.
        upperfencesrc
            Sets the source reference on Chart Studio Cloud for
            `upperfence`.
        visible
            Determines whether or not this trace is visible. If
            "legendonly", the trace is not drawn, but can appear as
            a legend item (provided that the legend itself is
            visible).
        whiskerwidth
            Sets the width of the whiskers relative to the box'
            width. For example, with 1, the whiskers are as wide as
            the box(es).
        width
            Sets the width of the box in data coordinate If 0
            (default value) the width is automatically selected
            based on the positions of other box traces in the same
            subplot.
        x
            Sets the x sample data or coordinates. See overview for
            more info.
        x0
            Sets the x coordinate for single-box traces or the
            starting coordinate for multi-box traces set using
            q1/median/q3. See overview for more info.
        xaxis
            Sets a reference between this trace's x coordinates and
            a 2D cartesian x axis. If "x" (the default value), the
            x coordinates refer to `layout.xaxis`. If "x2", the x
            coordinates refer to `layout.xaxis2`, and so on.
        xcalendar
            Sets the calendar system to use with `x` date data.
        xhoverformat
            Sets the hover text formatting rulefor `x`  using d3
            formatting mini-languages which are very similar to
            those in Python. For numbers, see:
            https://github.com/d3/d3-format/tree/v1.4.5#d3-format.
            And for dates see: https://github.com/d3/d3-time-
            format/tree/v2.2.3#locale_format. We add two items to
            d3's date formatter: "%h" for half of the year as a
            decimal number as well as "%{n}f" for fractional
            seconds with n digits. For example, *2016-10-13
            09:15:23.456* with tickformat "%H~%M~%S.%2f" would
            display *09~15~23.46*By default the values are
            formatted using `xaxis.hoverformat`.
        xperiod
            Only relevant when the axis `type` is "date". Sets the
            period positioning in milliseconds or "M<n>" on the x
            axis. Special values in the form of "M<n>" could be
            used to declare the number of months. In this case `n`
            must be a positive integer.
        xperiod0
            Only relevant when the axis `type` is "date". Sets the
            base for period positioning in milliseconds or date
            string on the x0 axis. When `x0period` is round number
            of weeks, the `x0period0` by default would be on a
            Sunday i.e. 2000-01-02, otherwise it would be at
            2000-01-01.
        xperiodalignment
            Only relevant when the axis `type` is "date". Sets the
            alignment of data points on the x axis.
        xsrc
            Sets the source reference on Chart Studio Cloud for
            `x`.
        y
            Sets the y sample data or coordinates. See overview for
            more info.
        y0
            Sets the y coordinate for single-box traces or the
            starting coordinate for multi-box traces set using
            q1/median/q3. See overview for more info.
        yaxis
            Sets a reference between this trace's y coordinates and
            a 2D cartesian y axis. If "y" (the default value), the
            y coordinates refer to `layout.yaxis`. If "y2", the y
            coordinates refer to `layout.yaxis2`, and so on.
        ycalendar
            Sets the calendar system to use with `y` date data.
        yhoverformat
            Sets the hover text formatting rulefor `y`  using d3
            formatting mini-languages which are very similar to
            those in Python. For numbers, see:
            https://github.com/d3/d3-format/tree/v1.4.5#d3-format.
            And for dates see: https://github.com/d3/d3-time-
            format/tree/v2.2.3#locale_format. We add two items to
            d3's date formatter: "%h" for half of the year as a
            decimal number as well as "%{n}f" for fractional
            seconds with n digits. For example, *2016-10-13
            09:15:23.456* with tickformat "%H~%M~%S.%2f" would
            display *09~15~23.46*By default the values are
            formatted using `yaxis.hoverformat`.
        yperiod
            Only relevant when the axis `type` is "date". Sets the
            period positioning in milliseconds or "M<n>" on the y
            axis. Special values in the form of "M<n>" could be
            used to declare the number of months. In this case `n`
            must be a positive integer.
        yperiod0
            Only relevant when the axis `type` is "date". Sets the
            base for period positioning in milliseconds or date
            string on the y0 axis. When `y0period` is round number
            of weeks, the `y0period0` by default would be on a
            Sunday i.e. 2000-01-02, otherwise it would be at
            2000-01-01.
        yperiodalignment
            Only relevant when the axis `type` is "date". Sets the
            alignment of data points on the y axis.
        ysrc
            Sets the source reference on Chart Studio Cloud for
            `y`.
        zorder
            Sets the layer on which this trace is displayed,
            relative to other SVG traces on the same subplot. SVG
            traces with higher `zorder` appear in front of those
            with lower `zorder`.
        """

    def __init__(
        self,
        arg=None,
        alignmentgroup=None,
        boxmean=None,
        boxpoints=None,
        customdata=None,
        customdatasrc=None,
        dx=None,
        dy=None,
        fillcolor=None,
        hoverinfo=None,
        hoverinfosrc=None,
        hoverlabel=None,
        hoveron=None,
        hovertemplate=None,
        hovertemplatesrc=None,
        hovertext=None,
        hovertextsrc=None,
        ids=None,
        idssrc=None,
        jitter=None,
        legend=None,
        legendgroup=None,
        legendgrouptitle=None,
        legendrank=None,
        legendwidth=None,
        line=None,
        lowerfence=None,
        lowerfencesrc=None,
        marker=None,
        mean=None,
        meansrc=None,
        median=None,
        mediansrc=None,
        meta=None,
        metasrc=None,
        name=None,
        notched=None,
        notchspan=None,
        notchspansrc=None,
        notchwidth=None,
        offsetgroup=None,
        opacity=None,
        orientation=None,
        pointpos=None,
        q1=None,
        q1src=None,
        q3=None,
        q3src=None,
        quartilemethod=None,
        sd=None,
        sdmultiple=None,
        sdsrc=None,
        selected=None,
        selectedpoints=None,
        showlegend=None,
        showwhiskers=None,
        sizemode=None,
        stream=None,
        text=None,
        textsrc=None,
        uid=None,
        uirevision=None,
        unselected=None,
        upperfence=None,
        upperfencesrc=None,
        visible=None,
        whiskerwidth=None,
        width=None,
        x=None,
        x0=None,
        xaxis=None,
        xcalendar=None,
        xhoverformat=None,
        xperiod=None,
        xperiod0=None,
        xperiodalignment=None,
        xsrc=None,
        y=None,
        y0=None,
        yaxis=None,
        ycalendar=None,
        yhoverformat=None,
        yperiod=None,
        yperiod0=None,
        yperiodalignment=None,
        ysrc=None,
        zorder=None,
        **kwargs,
    ):
        """
        Construct a new Box object

        Each box spans from quartile 1 (Q1) to quartile 3 (Q3). The
        second quartile (Q2, i.e. the median) is marked by a line
        inside the box. The fences grow outward from the boxes' edges,
        by default they span +/- 1.5 times the interquartile range
        (IQR: Q3-Q1), The sample mean and standard deviation as well as
        notches and the sample, outlier and suspected outliers points
        can be optionally added to the box plot. The values and
        positions corresponding to each boxes can be input using two
        signatures. The first signature expects users to supply the
        sample values in the `y` data array for vertical boxes (`x` for
        horizontal boxes). By supplying an `x` (`y`) array, one box per
        distinct `x` (`y`) value is drawn If no `x` (`y`) list is
        provided, a single box is drawn. In this case, the box is
        positioned with the trace `name` or with `x0` (`y0`) if
        provided. The second signature expects users to supply the
        boxes corresponding Q1, median and Q3 statistics in the `q1`,
        `median` and `q3` data arrays respectively. Other box features
        relying on statistics namely `lowerfence`, `upperfence`,
        `notchspan` can be set directly by the users. To have plotly
        compute them or to show sample points besides the boxes, users
        can set the `y` data array for vertical boxes (`x` for
        horizontal boxes) to a 2D array with the outer length
        corresponding to the number of boxes in the traces and the
        inner length corresponding the sample size.

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of :class:`plotly.graph_objs.Box`
        alignmentgroup
            Set several traces linked to the same position axis or
            matching axes to the same alignmentgroup. This controls
            whether bars compute their positional range dependently
            or independently.
        boxmean
            If True, the mean of the box(es)' underlying
            distribution is drawn as a dashed line inside the
            box(es). If "sd" the standard deviation is also drawn.
            Defaults to True when `mean` is set. Defaults to "sd"
            when `sd` is set Otherwise defaults to False.
        boxpoints
            If "outliers", only the sample points lying outside the
            whiskers are shown If "suspectedoutliers", the outlier
            points are shown and points either less than 4*Q1-3*Q3
            or greater than 4*Q3-3*Q1 are highlighted (see
            `outliercolor`) If "all", all sample points are shown
            If False, only the box(es) are shown with no sample
            points Defaults to "suspectedoutliers" when
            `marker.outliercolor` or `marker.line.outliercolor` is
            set. Defaults to "all" under the q1/median/q3
            signature. Otherwise defaults to "outliers".
        customdata
            Assigns extra data each datum. This may be useful when
            listening to hover, click and selection events. Note
            that, "scatter" traces also appends customdata items in
            the markers DOM elements
        customdatasrc
            Sets the source reference on Chart Studio Cloud for
            `customdata`.
        dx
            Sets the x coordinate step for multi-box traces set
            using q1/median/q3.
        dy
            Sets the y coordinate step for multi-box traces set
            using q1/median/q3.
        fillcolor
            Sets the fill color. Defaults to a half-transparent
            variant of the line color, marker color, or marker line
            color, whichever is available.
        hoverinfo
            Determines which trace information appear on hover. If
            `none` or `skip` are set, no information is displayed
            upon hovering. But, if `none` is set, click and hover
            events are still fired.
        hoverinfosrc
            Sets the source reference on Chart Studio Cloud for
            `hoverinfo`.
        hoverlabel
            :class:`plotly.graph_objects.box.Hoverlabel` instance
            or dict with compatible properties
        hoveron
            Do the hover effects highlight individual boxes  or
            sample points or both?
        hovertemplate
            Template string used for rendering the information that
            appear on hover box. Note that this will override
            `hoverinfo`. Variables are inserted using %{variable},
            for example "y: %{y}" as well as %{xother}, {%_xother},
            {%_xother_}, {%xother_}. When showing info for several
            points, "xother" will be added to those with different
            x positions from the first point. An underscore before
            or after "(x|y)other" will add a space on that side,
            only when this field is shown. Numbers are formatted
            using d3-format's syntax %{variable:d3-format}, for
            example "Price: %{y:$.2f}".
            https://github.com/d3/d3-format/tree/v1.4.5#d3-format
            for details on the formatting syntax. Dates are
            formatted using d3-time-format's syntax
            %{variable|d3-time-format}, for example "Day:
            %{2019-01-01|%A}". https://github.com/d3/d3-time-
            format/tree/v2.2.3#locale_format for details on the
            date formatting syntax. The variables available in
            `hovertemplate` are the ones emitted as event data
            described at this link
            https://plotly.com/javascript/plotlyjs-events/#event-
            data. Additionally, every attributes that can be
            specified per-point (the ones that are `arrayOk: true`)
            are available.  Anything contained in tag `<extra>` is
            displayed in the secondary box, for example
            `<extra>%{fullData.name}</extra>`. To hide the
            secondary box completely, use an empty tag
            `<extra></extra>`.
        hovertemplatesrc
            Sets the source reference on Chart Studio Cloud for
            `hovertemplate`.
        hovertext
            Same as `text`.
        hovertextsrc
            Sets the source reference on Chart Studio Cloud for
            `hovertext`.
        ids
            Assigns id labels to each datum. These ids for object
            constancy of data points during animation. Should be an
            array of strings, not numbers or any other type.
        idssrc
            Sets the source reference on Chart Studio Cloud for
            `ids`.
        jitter
            Sets the amount of jitter in the sample points drawn.
            If 0, the sample points align along the distribution
            axis. If 1, the sample points are drawn in a random
            jitter of width equal to the width of the box(es).
        legend
            Sets the reference to a legend to show this trace in.
            References to these legends are "legend", "legend2",
            "legend3", etc. Settings for these legends are set in
            the layout, under `layout.legend`, `layout.legend2`,
            etc.
        legendgroup
            Sets the legend group for this trace. Traces and shapes
            part of the same legend group hide/show at the same
            time when toggling legend items.
        legendgrouptitle
            :class:`plotly.graph_objects.box.Legendgrouptitle`
            instance or dict with compatible properties
        legendrank
            Sets the legend rank for this trace. Items and groups
            with smaller ranks are presented on top/left side while
            with "reversed" `legend.traceorder` they are on
            bottom/right side. The default legendrank is 1000, so
            that you can use ranks less than 1000 to place certain
            items before all unranked items, and ranks greater than
            1000 to go after all unranked items. When having
            unranked or equal rank items shapes would be displayed
            after traces i.e. according to their order in data and
            layout.
        legendwidth
            Sets the width (in px or fraction) of the legend for
            this trace.
        line
            :class:`plotly.graph_objects.box.Line` instance or dict
            with compatible properties
        lowerfence
            Sets the lower fence values. There should be as many
            items as the number of boxes desired. This attribute
            has effect only under the q1/median/q3 signature. If
            `lowerfence` is not provided but a sample (in `y` or
            `x`) is set, we compute the lower as the last sample
            point below 1.5 times the IQR.
        lowerfencesrc
            Sets the source reference on Chart Studio Cloud for
            `lowerfence`.
        marker
            :class:`plotly.graph_objects.box.Marker` instance or
            dict with compatible properties
        mean
            Sets the mean values. There should be as many items as
            the number of boxes desired. This attribute has effect
            only under the q1/median/q3 signature. If `mean` is not
            provided but a sample (in `y` or `x`) is set, we
            compute the mean for each box using the sample values.
        meansrc
            Sets the source reference on Chart Studio Cloud for
            `mean`.
        median
            Sets the median values. There should be as many items
            as the number of boxes desired.
        mediansrc
            Sets the source reference on Chart Studio Cloud for
            `median`.
        meta
            Assigns extra meta information associated with this
            trace that can be used in various text attributes.
            Attributes such as trace `name`, graph, axis and
            colorbar `title.text`, annotation `text`
            `rangeselector`, `updatemenues` and `sliders` `label`
            text all support `meta`. To access the trace `meta`
            values in an attribute in the same trace, simply use
            `%{meta[i]}` where `i` is the index or key of the
            `meta` item in question. To access trace `meta` in
            layout attributes, use `%{data[n[.meta[i]}` where `i`
            is the index or key of the `meta` and `n` is the trace
            index.
        metasrc
            Sets the source reference on Chart Studio Cloud for
            `meta`.
        name
            Sets the trace name. The trace name appears as the
            legend item and on hover. For box traces, the name will
            also be used for the position coordinate, if `x` and
            `x0` (`y` and `y0` if horizontal) are missing and the
            position axis is categorical
        notched
            Determines whether or not notches are drawn. Notches
            displays a confidence interval around the median. We
            compute the confidence interval as median +/- 1.57 *
            IQR / sqrt(N), where IQR is the interquartile range and
            N is the sample size. If two boxes' notches do not
            overlap there is 95% confidence their medians differ.
            See https://sites.google.com/site/davidsstatistics/home
            /notched-box-plots for more info. Defaults to False
            unless `notchwidth` or `notchspan` is set.
        notchspan
            Sets the notch span from the boxes' `median` values.
            There should be as many items as the number of boxes
            desired. This attribute has effect only under the
            q1/median/q3 signature. If `notchspan` is not provided
            but a sample (in `y` or `x`) is set, we compute it as
            1.57 * IQR / sqrt(N), where N is the sample size.
        notchspansrc
            Sets the source reference on Chart Studio Cloud for
            `notchspan`.
        notchwidth
            Sets the width of the notches relative to the box'
            width. For example, with 0, the notches are as wide as
            the box(es).
        offsetgroup
            Set several traces linked to the same position axis or
            matching axes to the same offsetgroup where bars of the
            same position coordinate will line up.
        opacity
            Sets the opacity of the trace.
        orientation
            Sets the orientation of the box(es). If "v" ("h"), the
            distribution is visualized along the vertical
            (horizontal).
        pointpos
            Sets the position of the sample points in relation to
            the box(es). If 0, the sample points are places over
            the center of the box(es). Positive (negative) values
            correspond to positions to the right (left) for
            vertical boxes and above (below) for horizontal boxes
        q1
            Sets the Quartile 1 values. There should be as many
            items as the number of boxes desired.
        q1src
            Sets the source reference on Chart Studio Cloud for
            `q1`.
        q3
            Sets the Quartile 3 values. There should be as many
            items as the number of boxes desired.
        q3src
            Sets the source reference on Chart Studio Cloud for
            `q3`.
        quartilemethod
            Sets the method used to compute the sample's Q1 and Q3
            quartiles. The "linear" method uses the 25th percentile
            for Q1 and 75th percentile for Q3 as computed using
            method #10 (listed on
            http://jse.amstat.org/v14n3/langford.html). The
            "exclusive" method uses the median to divide the
            ordered dataset into two halves if the sample is odd,
            it does not include the median in either half - Q1 is
            then the median of the lower half and Q3 the median of
            the upper half. The "inclusive" method also uses the
            median to divide the ordered dataset into two halves
            but if the sample is odd, it includes the median in
            both halves - Q1 is then the median of the lower half
            and Q3 the median of the upper half.
        sd
            Sets the standard deviation values. There should be as
            many items as the number of boxes desired. This
            attribute has effect only under the q1/median/q3
            signature. If `sd` is not provided but a sample (in `y`
            or `x`) is set, we compute the standard deviation for
            each box using the sample values.
        sdmultiple
            Scales the box size when sizemode=sd Allowing boxes to
            be drawn across any stddev range For example 1-stddev,
            3-stddev, 5-stddev
        sdsrc
            Sets the source reference on Chart Studio Cloud for
            `sd`.
        selected
            :class:`plotly.graph_objects.box.Selected` instance or
            dict with compatible properties
        selectedpoints
            Array containing integer indices of selected points.
            Has an effect only for traces that support selections.
            Note that an empty array means an empty selection where
            the `unselected` are turned on for all points, whereas,
            any other non-array values means no selection all where
            the `selected` and `unselected` styles have no effect.
        showlegend
            Determines whether or not an item corresponding to this
            trace is shown in the legend.
        showwhiskers
            Determines whether or not whiskers are visible.
            Defaults to true for `sizemode` "quartiles", false for
            "sd".
        sizemode
            Sets the upper and lower bound for the boxes quartiles
            means box is drawn between Q1 and Q3 SD means the box
            is drawn between Mean +- Standard Deviation Argument
            sdmultiple (default 1) to scale the box size So it
            could be drawn 1-stddev, 3-stddev etc
        stream
            :class:`plotly.graph_objects.box.Stream` instance or
            dict with compatible properties
        text
            Sets the text elements associated with each sample
            value. If a single string, the same string appears over
            all the data points. If an array of string, the items
            are mapped in order to the this trace's (x,y)
            coordinates. To be seen, trace `hoverinfo` must contain
            a "text" flag.
        textsrc
            Sets the source reference on Chart Studio Cloud for
            `text`.
        uid
            Assign an id to this trace, Use this to provide object
            constancy between traces during animations and
            transitions.
        uirevision
            Controls persistence of some user-driven changes to the
            trace: `constraintrange` in `parcoords` traces, as well
            as some `editable: true` modifications such as `name`
            and `colorbar.title`. Defaults to `layout.uirevision`.
            Note that other user-driven trace attribute changes are
            controlled by `layout` attributes: `trace.visible` is
            controlled by `layout.legend.uirevision`,
            `selectedpoints` is controlled by
            `layout.selectionrevision`, and `colorbar.(x|y)`
            (accessible with `config: {editable: true}`) is
            controlled by `layout.editrevision`. Trace changes are
            tracked by `uid`, which only falls back on trace index
            if no `uid` is provided. So if your app can add/remove
            traces before the end of the `data` array, such that
            the same trace has a different index, you can still
            preserve user-driven changes if you give each trace a
            `uid` that stays with it as it moves.
        unselected
            :class:`plotly.graph_objects.box.Unselected` instance
            or dict with compatible properties
        upperfence
            Sets the upper fence values. There should be as many
            items as the number of boxes desired. This attribute
            has effect only under the q1/median/q3 signature. If
            `upperfence` is not provided but a sample (in `y` or
            `x`) is set, we compute the upper as the last sample
            point above 1.5 times the IQR.
        upperfencesrc
            Sets the source reference on Chart Studio Cloud for
            `upperfence`.
        visible
            Determines whether or not this trace is visible. If
            "legendonly", the trace is not drawn, but can appear as
            a legend item (provided that the legend itself is
            visible).
        whiskerwidth
            Sets the width of the whiskers relative to the box'
            width. For example, with 1, the whiskers are as wide as
            the box(es).
        width
            Sets the width of the box in data coordinate If 0
            (default value) the width is automatically selected
            based on the positions of other box traces in the same
            subplot.
        x
            Sets the x sample data or coordinates. See overview for
            more info.
        x0
            Sets the x coordinate for single-box traces or the
            starting coordinate for multi-box traces set using
            q1/median/q3. See overview for more info.
        xaxis
            Sets a reference between this trace's x coordinates and
            a 2D cartesian x axis. If "x" (the default value), the
            x coordinates refer to `layout.xaxis`. If "x2", the x
            coordinates refer to `layout.xaxis2`, and so on.
        xcalendar
            Sets the calendar system to use with `x` date data.
        xhoverformat
            Sets the hover text formatting rulefor `x`  using d3
            formatting mini-languages which are very similar to
            those in Python. For numbers, see:
            https://github.com/d3/d3-format/tree/v1.4.5#d3-format.
            And for dates see: https://github.com/d3/d3-time-
            format/tree/v2.2.3#locale_format. We add two items to
            d3's date formatter: "%h" for half of the year as a
            decimal number as well as "%{n}f" for fractional
            seconds with n digits. For example, *2016-10-13
            09:15:23.456* with tickformat "%H~%M~%S.%2f" would
            display *09~15~23.46*By default the values are
            formatted using `xaxis.hoverformat`.
        xperiod
            Only relevant when the axis `type` is "date". Sets the
            period positioning in milliseconds or "M<n>" on the x
            axis. Special values in the form of "M<n>" could be
            used to declare the number of months. In this case `n`
            must be a positive integer.
        xperiod0
            Only relevant when the axis `type` is "date". Sets the
            base for period positioning in milliseconds or date
            string on the x0 axis. When `x0period` is round number
            of weeks, the `x0period0` by default would be on a
            Sunday i.e. 2000-01-02, otherwise it would be at
            2000-01-01.
        xperiodalignment
            Only relevant when the axis `type` is "date". Sets the
            alignment of data points on the x axis.
        xsrc
            Sets the source reference on Chart Studio Cloud for
            `x`.
        y
            Sets the y sample data or coordinates. See overview for
            more info.
        y0
            Sets the y coordinate for single-box traces or the
            starting coordinate for multi-box traces set using
            q1/median/q3. See overview for more info.
        yaxis
            Sets a reference between this trace's y coordinates and
            a 2D cartesian y axis. If "y" (the default value), the
            y coordinates refer to `layout.yaxis`. If "y2", the y
            coordinates refer to `layout.yaxis2`, and so on.
        ycalendar
            Sets the calendar system to use with `y` date data.
        yhoverformat
            Sets the hover text formatting rulefor `y`  using d3
            formatting mini-languages which are very similar to
            those in Python. For numbers, see:
            https://github.com/d3/d3-format/tree/v1.4.5#d3-format.
            And for dates see: https://github.com/d3/d3-time-
            format/tree/v2.2.3#locale_format. We add two items to
            d3's date formatter: "%h" for half of the year as a
            decimal number as well as "%{n}f" for fractional
            seconds with n digits. For example, *2016-10-13
            09:15:23.456* with tickformat "%H~%M~%S.%2f" would
            display *09~15~23.46*By default the values are
            formatted using `yaxis.hoverformat`.
        yperiod
            Only relevant when the axis `type` is "date". Sets the
            period positioning in milliseconds or "M<n>" on the y
            axis. Special values in the form of "M<n>" could be
            used to declare the number of months. In this case `n`
            must be a positive integer.
        yperiod0
            Only relevant when the axis `type` is "date". Sets the
            base for period positioning in milliseconds or date
            string on the y0 axis. When `y0period` is round number
            of weeks, the `y0period0` by default would be on a
            Sunday i.e. 2000-01-02, otherwise it would be at
            2000-01-01.
        yperiodalignment
            Only relevant when the axis `type` is "date". Sets the
            alignment of data points on the y axis.
        ysrc
            Sets the source reference on Chart Studio Cloud for
            `y`.
        zorder
            Sets the layer on which this trace is displayed,
            relative to other SVG traces on the same subplot. SVG
            traces with higher `zorder` appear in front of those
            with lower `zorder`.

        Returns
        -------
        Box
        """
        super().__init__("box")
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
The first argument to the plotly.graph_objs.Box
constructor must be a dict or
an instance of :class:`plotly.graph_objs.Box`""")

        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        self._set_property("alignmentgroup", arg, alignmentgroup)
        self._set_property("boxmean", arg, boxmean)
        self._set_property("boxpoints", arg, boxpoints)
        self._set_property("customdata", arg, customdata)
        self._set_property("customdatasrc", arg, customdatasrc)
        self._set_property("dx", arg, dx)
        self._set_property("dy", arg, dy)
        self._set_property("fillcolor", arg, fillcolor)
        self._set_property("hoverinfo", arg, hoverinfo)
        self._set_property("hoverinfosrc", arg, hoverinfosrc)
        self._set_property("hoverlabel", arg, hoverlabel)
        self._set_property("hoveron", arg, hoveron)
        self._set_property("hovertemplate", arg, hovertemplate)
        self._set_property("hovertemplatesrc", arg, hovertemplatesrc)
        self._set_property("hovertext", arg, hovertext)
        self._set_property("hovertextsrc", arg, hovertextsrc)
        self._set_property("ids", arg, ids)
        self._set_property("idssrc", arg, idssrc)
        self._set_property("jitter", arg, jitter)
        self._set_property("legend", arg, legend)
        self._set_property("legendgroup", arg, legendgroup)
        self._set_property("legendgrouptitle", arg, legendgrouptitle)
        self._set_property("legendrank", arg, legendrank)
        self._set_property("legendwidth", arg, legendwidth)
        self._set_property("line", arg, line)
        self._set_property("lowerfence", arg, lowerfence)
        self._set_property("lowerfencesrc", arg, lowerfencesrc)
        self._set_property("marker", arg, marker)
        self._set_property("mean", arg, mean)
        self._set_property("meansrc", arg, meansrc)
        self._set_property("median", arg, median)
        self._set_property("mediansrc", arg, mediansrc)
        self._set_property("meta", arg, meta)
        self._set_property("metasrc", arg, metasrc)
        self._set_property("name", arg, name)
        self._set_property("notched", arg, notched)
        self._set_property("notchspan", arg, notchspan)
        self._set_property("notchspansrc", arg, notchspansrc)
        self._set_property("notchwidth", arg, notchwidth)
        self._set_property("offsetgroup", arg, offsetgroup)
        self._set_property("opacity", arg, opacity)
        self._set_property("orientation", arg, orientation)
        self._set_property("pointpos", arg, pointpos)
        self._set_property("q1", arg, q1)
        self._set_property("q1src", arg, q1src)
        self._set_property("q3", arg, q3)
        self._set_property("q3src", arg, q3src)
        self._set_property("quartilemethod", arg, quartilemethod)
        self._set_property("sd", arg, sd)
        self._set_property("sdmultiple", arg, sdmultiple)
        self._set_property("sdsrc", arg, sdsrc)
        self._set_property("selected", arg, selected)
        self._set_property("selectedpoints", arg, selectedpoints)
        self._set_property("showlegend", arg, showlegend)
        self._set_property("showwhiskers", arg, showwhiskers)
        self._set_property("sizemode", arg, sizemode)
        self._set_property("stream", arg, stream)
        self._set_property("text", arg, text)
        self._set_property("textsrc", arg, textsrc)
        self._set_property("uid", arg, uid)
        self._set_property("uirevision", arg, uirevision)
        self._set_property("unselected", arg, unselected)
        self._set_property("upperfence", arg, upperfence)
        self._set_property("upperfencesrc", arg, upperfencesrc)
        self._set_property("visible", arg, visible)
        self._set_property("whiskerwidth", arg, whiskerwidth)
        self._set_property("width", arg, width)
        self._set_property("x", arg, x)
        self._set_property("x0", arg, x0)
        self._set_property("xaxis", arg, xaxis)
        self._set_property("xcalendar", arg, xcalendar)
        self._set_property("xhoverformat", arg, xhoverformat)
        self._set_property("xperiod", arg, xperiod)
        self._set_property("xperiod0", arg, xperiod0)
        self._set_property("xperiodalignment", arg, xperiodalignment)
        self._set_property("xsrc", arg, xsrc)
        self._set_property("y", arg, y)
        self._set_property("y0", arg, y0)
        self._set_property("yaxis", arg, yaxis)
        self._set_property("ycalendar", arg, ycalendar)
        self._set_property("yhoverformat", arg, yhoverformat)
        self._set_property("yperiod", arg, yperiod)
        self._set_property("yperiod0", arg, yperiod0)
        self._set_property("yperiodalignment", arg, yperiodalignment)
        self._set_property("ysrc", arg, ysrc)
        self._set_property("zorder", arg, zorder)

        self._props["type"] = "box"
        arg.pop("type", None)
        self._process_kwargs(**dict(arg, **kwargs))
        self._skip_invalid = False
