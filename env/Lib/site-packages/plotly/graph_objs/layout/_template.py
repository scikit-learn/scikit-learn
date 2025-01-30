from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Template(_BaseLayoutHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "layout"
    _path_str = "layout.template"
    _valid_props = {"data", "layout"}

    # data
    # ----
    @property
    def data(self):
        """
        The 'data' property is an instance of Data
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.layout.template.Data`
          - A dict of string/value properties that will be passed
            to the Data constructor

            Supported dict properties:

                barpolar
                    A tuple of
                    :class:`plotly.graph_objects.Barpolar`
                    instances or dicts with compatible properties
                bar
                    A tuple of :class:`plotly.graph_objects.Bar`
                    instances or dicts with compatible properties
                box
                    A tuple of :class:`plotly.graph_objects.Box`
                    instances or dicts with compatible properties
                candlestick
                    A tuple of
                    :class:`plotly.graph_objects.Candlestick`
                    instances or dicts with compatible properties
                carpet
                    A tuple of :class:`plotly.graph_objects.Carpet`
                    instances or dicts with compatible properties
                choroplethmapbox
                    A tuple of
                    :class:`plotly.graph_objects.Choroplethmapbox`
                    instances or dicts with compatible properties
                choroplethmap
                    A tuple of
                    :class:`plotly.graph_objects.Choroplethmap`
                    instances or dicts with compatible properties
                choropleth
                    A tuple of
                    :class:`plotly.graph_objects.Choropleth`
                    instances or dicts with compatible properties
                cone
                    A tuple of :class:`plotly.graph_objects.Cone`
                    instances or dicts with compatible properties
                contourcarpet
                    A tuple of
                    :class:`plotly.graph_objects.Contourcarpet`
                    instances or dicts with compatible properties
                contour
                    A tuple of
                    :class:`plotly.graph_objects.Contour` instances
                    or dicts with compatible properties
                densitymapbox
                    A tuple of
                    :class:`plotly.graph_objects.Densitymapbox`
                    instances or dicts with compatible properties
                densitymap
                    A tuple of
                    :class:`plotly.graph_objects.Densitymap`
                    instances or dicts with compatible properties
                funnelarea
                    A tuple of
                    :class:`plotly.graph_objects.Funnelarea`
                    instances or dicts with compatible properties
                funnel
                    A tuple of :class:`plotly.graph_objects.Funnel`
                    instances or dicts with compatible properties
                heatmapgl
                    A tuple of
                    :class:`plotly.graph_objects.Heatmapgl`
                    instances or dicts with compatible properties
                heatmap
                    A tuple of
                    :class:`plotly.graph_objects.Heatmap` instances
                    or dicts with compatible properties
                histogram2dcontour
                    A tuple of :class:`plotly.graph_objects.Histogr
                    am2dContour` instances or dicts with compatible
                    properties
                histogram2d
                    A tuple of
                    :class:`plotly.graph_objects.Histogram2d`
                    instances or dicts with compatible properties
                histogram
                    A tuple of
                    :class:`plotly.graph_objects.Histogram`
                    instances or dicts with compatible properties
                icicle
                    A tuple of :class:`plotly.graph_objects.Icicle`
                    instances or dicts with compatible properties
                image
                    A tuple of :class:`plotly.graph_objects.Image`
                    instances or dicts with compatible properties
                indicator
                    A tuple of
                    :class:`plotly.graph_objects.Indicator`
                    instances or dicts with compatible properties
                isosurface
                    A tuple of
                    :class:`plotly.graph_objects.Isosurface`
                    instances or dicts with compatible properties
                mesh3d
                    A tuple of :class:`plotly.graph_objects.Mesh3d`
                    instances or dicts with compatible properties
                ohlc
                    A tuple of :class:`plotly.graph_objects.Ohlc`
                    instances or dicts with compatible properties
                parcats
                    A tuple of
                    :class:`plotly.graph_objects.Parcats` instances
                    or dicts with compatible properties
                parcoords
                    A tuple of
                    :class:`plotly.graph_objects.Parcoords`
                    instances or dicts with compatible properties
                pie
                    A tuple of :class:`plotly.graph_objects.Pie`
                    instances or dicts with compatible properties
                pointcloud
                    A tuple of
                    :class:`plotly.graph_objects.Pointcloud`
                    instances or dicts with compatible properties
                sankey
                    A tuple of :class:`plotly.graph_objects.Sankey`
                    instances or dicts with compatible properties
                scatter3d
                    A tuple of
                    :class:`plotly.graph_objects.Scatter3d`
                    instances or dicts with compatible properties
                scattercarpet
                    A tuple of
                    :class:`plotly.graph_objects.Scattercarpet`
                    instances or dicts with compatible properties
                scattergeo
                    A tuple of
                    :class:`plotly.graph_objects.Scattergeo`
                    instances or dicts with compatible properties
                scattergl
                    A tuple of
                    :class:`plotly.graph_objects.Scattergl`
                    instances or dicts with compatible properties
                scattermapbox
                    A tuple of
                    :class:`plotly.graph_objects.Scattermapbox`
                    instances or dicts with compatible properties
                scattermap
                    A tuple of
                    :class:`plotly.graph_objects.Scattermap`
                    instances or dicts with compatible properties
                scatterpolargl
                    A tuple of
                    :class:`plotly.graph_objects.Scatterpolargl`
                    instances or dicts with compatible properties
                scatterpolar
                    A tuple of
                    :class:`plotly.graph_objects.Scatterpolar`
                    instances or dicts with compatible properties
                scatter
                    A tuple of
                    :class:`plotly.graph_objects.Scatter` instances
                    or dicts with compatible properties
                scattersmith
                    A tuple of
                    :class:`plotly.graph_objects.Scattersmith`
                    instances or dicts with compatible properties
                scatterternary
                    A tuple of
                    :class:`plotly.graph_objects.Scatterternary`
                    instances or dicts with compatible properties
                splom
                    A tuple of :class:`plotly.graph_objects.Splom`
                    instances or dicts with compatible properties
                streamtube
                    A tuple of
                    :class:`plotly.graph_objects.Streamtube`
                    instances or dicts with compatible properties
                sunburst
                    A tuple of
                    :class:`plotly.graph_objects.Sunburst`
                    instances or dicts with compatible properties
                surface
                    A tuple of
                    :class:`plotly.graph_objects.Surface` instances
                    or dicts with compatible properties
                table
                    A tuple of :class:`plotly.graph_objects.Table`
                    instances or dicts with compatible properties
                treemap
                    A tuple of
                    :class:`plotly.graph_objects.Treemap` instances
                    or dicts with compatible properties
                violin
                    A tuple of :class:`plotly.graph_objects.Violin`
                    instances or dicts with compatible properties
                volume
                    A tuple of :class:`plotly.graph_objects.Volume`
                    instances or dicts with compatible properties
                waterfall
                    A tuple of
                    :class:`plotly.graph_objects.Waterfall`
                    instances or dicts with compatible properties

        Returns
        -------
        plotly.graph_objs.layout.template.Data
        """
        return self["data"]

    @data.setter
    def data(self, val):
        self["data"] = val

    # layout
    # ------
    @property
    def layout(self):
        """
        The 'layout' property is an instance of Layout
        that may be specified as:
          - An instance of :class:`plotly.graph_objs.Layout`
          - A dict of string/value properties that will be passed
            to the Layout constructor

            Supported dict properties:

        Returns
        -------
        plotly.graph_objs.layout.template.Layout
        """
        return self["layout"]

    @layout.setter
    def layout(self, val):
        self["layout"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        data
            :class:`plotly.graph_objects.layout.template.Data`
            instance or dict with compatible properties
        layout
            :class:`plotly.graph_objects.Layout` instance or dict
            with compatible properties
        """

    def __init__(self, arg=None, data=None, layout=None, **kwargs):
        """
        Construct a new Template object

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

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.layout.Template`
        data
            :class:`plotly.graph_objects.layout.template.Data`
            instance or dict with compatible properties
        layout
            :class:`plotly.graph_objects.Layout` instance or dict
            with compatible properties

        Returns
        -------
        Template
        """
        super(Template, self).__init__("template")

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
The first argument to the plotly.graph_objs.layout.Template
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.Template`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("data", None)
        _v = data if data is not None else _v
        if _v is not None:
            self["data"] = _v
        _v = arg.pop("layout", None)
        _v = layout if layout is not None else _v
        if _v is not None:
            self["layout"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
