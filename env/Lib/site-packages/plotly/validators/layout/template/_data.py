import _plotly_utils.basevalidators


class DataValidator(_plotly_utils.basevalidators.CompoundValidator):
    def __init__(self, plotly_name="data", parent_name="layout.template", **kwargs):
        super(DataValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=kwargs.pop("data_class_str", "Data"),
            data_docs=kwargs.pop(
                "data_docs",
                """
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
""",
            ),
            **kwargs,
        )
