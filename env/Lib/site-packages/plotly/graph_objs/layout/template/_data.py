from plotly.basedatatypes import BaseLayoutHierarchyType as _BaseLayoutHierarchyType
import copy as _copy


class Data(_BaseLayoutHierarchyType):

    # class properties
    # --------------------
    _parent_path_str = "layout.template"
    _path_str = "layout.template.data"
    _valid_props = {
        "bar",
        "barpolar",
        "box",
        "candlestick",
        "carpet",
        "choropleth",
        "choroplethmap",
        "choroplethmapbox",
        "cone",
        "contour",
        "contourcarpet",
        "densitymap",
        "densitymapbox",
        "funnel",
        "funnelarea",
        "heatmap",
        "heatmapgl",
        "histogram",
        "histogram2d",
        "histogram2dcontour",
        "icicle",
        "image",
        "indicator",
        "isosurface",
        "mesh3d",
        "ohlc",
        "parcats",
        "parcoords",
        "pie",
        "pointcloud",
        "sankey",
        "scatter",
        "scatter3d",
        "scattercarpet",
        "scattergeo",
        "scattergl",
        "scattermap",
        "scattermapbox",
        "scatterpolar",
        "scatterpolargl",
        "scattersmith",
        "scatterternary",
        "splom",
        "streamtube",
        "sunburst",
        "surface",
        "table",
        "treemap",
        "violin",
        "volume",
        "waterfall",
    }

    # barpolar
    # --------
    @property
    def barpolar(self):
        """
        The 'barpolar' property is a tuple of instances of
        Barpolar that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Barpolar
          - A list or tuple of dicts of string/value properties that
            will be passed to the Barpolar constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Barpolar]
        """
        return self["barpolar"]

    @barpolar.setter
    def barpolar(self, val):
        self["barpolar"] = val

    # bar
    # ---
    @property
    def bar(self):
        """
        The 'bar' property is a tuple of instances of
        Bar that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Bar
          - A list or tuple of dicts of string/value properties that
            will be passed to the Bar constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Bar]
        """
        return self["bar"]

    @bar.setter
    def bar(self, val):
        self["bar"] = val

    # box
    # ---
    @property
    def box(self):
        """
        The 'box' property is a tuple of instances of
        Box that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Box
          - A list or tuple of dicts of string/value properties that
            will be passed to the Box constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Box]
        """
        return self["box"]

    @box.setter
    def box(self, val):
        self["box"] = val

    # candlestick
    # -----------
    @property
    def candlestick(self):
        """
        The 'candlestick' property is a tuple of instances of
        Candlestick that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Candlestick
          - A list or tuple of dicts of string/value properties that
            will be passed to the Candlestick constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Candlestick]
        """
        return self["candlestick"]

    @candlestick.setter
    def candlestick(self, val):
        self["candlestick"] = val

    # carpet
    # ------
    @property
    def carpet(self):
        """
        The 'carpet' property is a tuple of instances of
        Carpet that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Carpet
          - A list or tuple of dicts of string/value properties that
            will be passed to the Carpet constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Carpet]
        """
        return self["carpet"]

    @carpet.setter
    def carpet(self, val):
        self["carpet"] = val

    # choroplethmapbox
    # ----------------
    @property
    def choroplethmapbox(self):
        """
        The 'choroplethmapbox' property is a tuple of instances of
        Choroplethmapbox that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Choroplethmapbox
          - A list or tuple of dicts of string/value properties that
            will be passed to the Choroplethmapbox constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Choroplethmapbox]
        """
        return self["choroplethmapbox"]

    @choroplethmapbox.setter
    def choroplethmapbox(self, val):
        self["choroplethmapbox"] = val

    # choroplethmap
    # -------------
    @property
    def choroplethmap(self):
        """
        The 'choroplethmap' property is a tuple of instances of
        Choroplethmap that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Choroplethmap
          - A list or tuple of dicts of string/value properties that
            will be passed to the Choroplethmap constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Choroplethmap]
        """
        return self["choroplethmap"]

    @choroplethmap.setter
    def choroplethmap(self, val):
        self["choroplethmap"] = val

    # choropleth
    # ----------
    @property
    def choropleth(self):
        """
        The 'choropleth' property is a tuple of instances of
        Choropleth that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Choropleth
          - A list or tuple of dicts of string/value properties that
            will be passed to the Choropleth constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Choropleth]
        """
        return self["choropleth"]

    @choropleth.setter
    def choropleth(self, val):
        self["choropleth"] = val

    # cone
    # ----
    @property
    def cone(self):
        """
        The 'cone' property is a tuple of instances of
        Cone that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Cone
          - A list or tuple of dicts of string/value properties that
            will be passed to the Cone constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Cone]
        """
        return self["cone"]

    @cone.setter
    def cone(self, val):
        self["cone"] = val

    # contourcarpet
    # -------------
    @property
    def contourcarpet(self):
        """
        The 'contourcarpet' property is a tuple of instances of
        Contourcarpet that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Contourcarpet
          - A list or tuple of dicts of string/value properties that
            will be passed to the Contourcarpet constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Contourcarpet]
        """
        return self["contourcarpet"]

    @contourcarpet.setter
    def contourcarpet(self, val):
        self["contourcarpet"] = val

    # contour
    # -------
    @property
    def contour(self):
        """
        The 'contour' property is a tuple of instances of
        Contour that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Contour
          - A list or tuple of dicts of string/value properties that
            will be passed to the Contour constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Contour]
        """
        return self["contour"]

    @contour.setter
    def contour(self, val):
        self["contour"] = val

    # densitymapbox
    # -------------
    @property
    def densitymapbox(self):
        """
        The 'densitymapbox' property is a tuple of instances of
        Densitymapbox that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Densitymapbox
          - A list or tuple of dicts of string/value properties that
            will be passed to the Densitymapbox constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Densitymapbox]
        """
        return self["densitymapbox"]

    @densitymapbox.setter
    def densitymapbox(self, val):
        self["densitymapbox"] = val

    # densitymap
    # ----------
    @property
    def densitymap(self):
        """
        The 'densitymap' property is a tuple of instances of
        Densitymap that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Densitymap
          - A list or tuple of dicts of string/value properties that
            will be passed to the Densitymap constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Densitymap]
        """
        return self["densitymap"]

    @densitymap.setter
    def densitymap(self, val):
        self["densitymap"] = val

    # funnelarea
    # ----------
    @property
    def funnelarea(self):
        """
        The 'funnelarea' property is a tuple of instances of
        Funnelarea that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Funnelarea
          - A list or tuple of dicts of string/value properties that
            will be passed to the Funnelarea constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Funnelarea]
        """
        return self["funnelarea"]

    @funnelarea.setter
    def funnelarea(self, val):
        self["funnelarea"] = val

    # funnel
    # ------
    @property
    def funnel(self):
        """
        The 'funnel' property is a tuple of instances of
        Funnel that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Funnel
          - A list or tuple of dicts of string/value properties that
            will be passed to the Funnel constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Funnel]
        """
        return self["funnel"]

    @funnel.setter
    def funnel(self, val):
        self["funnel"] = val

    # heatmapgl
    # ---------
    @property
    def heatmapgl(self):
        """
        The 'heatmapgl' property is a tuple of instances of
        Heatmapgl that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Heatmapgl
          - A list or tuple of dicts of string/value properties that
            will be passed to the Heatmapgl constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Heatmapgl]
        """
        return self["heatmapgl"]

    @heatmapgl.setter
    def heatmapgl(self, val):
        self["heatmapgl"] = val

    # heatmap
    # -------
    @property
    def heatmap(self):
        """
        The 'heatmap' property is a tuple of instances of
        Heatmap that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Heatmap
          - A list or tuple of dicts of string/value properties that
            will be passed to the Heatmap constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Heatmap]
        """
        return self["heatmap"]

    @heatmap.setter
    def heatmap(self, val):
        self["heatmap"] = val

    # histogram2dcontour
    # ------------------
    @property
    def histogram2dcontour(self):
        """
        The 'histogram2dcontour' property is a tuple of instances of
        Histogram2dContour that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Histogram2dContour
          - A list or tuple of dicts of string/value properties that
            will be passed to the Histogram2dContour constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Histogram2dContour]
        """
        return self["histogram2dcontour"]

    @histogram2dcontour.setter
    def histogram2dcontour(self, val):
        self["histogram2dcontour"] = val

    # histogram2d
    # -----------
    @property
    def histogram2d(self):
        """
        The 'histogram2d' property is a tuple of instances of
        Histogram2d that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Histogram2d
          - A list or tuple of dicts of string/value properties that
            will be passed to the Histogram2d constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Histogram2d]
        """
        return self["histogram2d"]

    @histogram2d.setter
    def histogram2d(self, val):
        self["histogram2d"] = val

    # histogram
    # ---------
    @property
    def histogram(self):
        """
        The 'histogram' property is a tuple of instances of
        Histogram that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Histogram
          - A list or tuple of dicts of string/value properties that
            will be passed to the Histogram constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Histogram]
        """
        return self["histogram"]

    @histogram.setter
    def histogram(self, val):
        self["histogram"] = val

    # icicle
    # ------
    @property
    def icicle(self):
        """
        The 'icicle' property is a tuple of instances of
        Icicle that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Icicle
          - A list or tuple of dicts of string/value properties that
            will be passed to the Icicle constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Icicle]
        """
        return self["icicle"]

    @icicle.setter
    def icicle(self, val):
        self["icicle"] = val

    # image
    # -----
    @property
    def image(self):
        """
        The 'image' property is a tuple of instances of
        Image that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Image
          - A list or tuple of dicts of string/value properties that
            will be passed to the Image constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Image]
        """
        return self["image"]

    @image.setter
    def image(self, val):
        self["image"] = val

    # indicator
    # ---------
    @property
    def indicator(self):
        """
        The 'indicator' property is a tuple of instances of
        Indicator that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Indicator
          - A list or tuple of dicts of string/value properties that
            will be passed to the Indicator constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Indicator]
        """
        return self["indicator"]

    @indicator.setter
    def indicator(self, val):
        self["indicator"] = val

    # isosurface
    # ----------
    @property
    def isosurface(self):
        """
        The 'isosurface' property is a tuple of instances of
        Isosurface that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Isosurface
          - A list or tuple of dicts of string/value properties that
            will be passed to the Isosurface constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Isosurface]
        """
        return self["isosurface"]

    @isosurface.setter
    def isosurface(self, val):
        self["isosurface"] = val

    # mesh3d
    # ------
    @property
    def mesh3d(self):
        """
        The 'mesh3d' property is a tuple of instances of
        Mesh3d that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Mesh3d
          - A list or tuple of dicts of string/value properties that
            will be passed to the Mesh3d constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Mesh3d]
        """
        return self["mesh3d"]

    @mesh3d.setter
    def mesh3d(self, val):
        self["mesh3d"] = val

    # ohlc
    # ----
    @property
    def ohlc(self):
        """
        The 'ohlc' property is a tuple of instances of
        Ohlc that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Ohlc
          - A list or tuple of dicts of string/value properties that
            will be passed to the Ohlc constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Ohlc]
        """
        return self["ohlc"]

    @ohlc.setter
    def ohlc(self, val):
        self["ohlc"] = val

    # parcats
    # -------
    @property
    def parcats(self):
        """
        The 'parcats' property is a tuple of instances of
        Parcats that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Parcats
          - A list or tuple of dicts of string/value properties that
            will be passed to the Parcats constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Parcats]
        """
        return self["parcats"]

    @parcats.setter
    def parcats(self, val):
        self["parcats"] = val

    # parcoords
    # ---------
    @property
    def parcoords(self):
        """
        The 'parcoords' property is a tuple of instances of
        Parcoords that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Parcoords
          - A list or tuple of dicts of string/value properties that
            will be passed to the Parcoords constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Parcoords]
        """
        return self["parcoords"]

    @parcoords.setter
    def parcoords(self, val):
        self["parcoords"] = val

    # pie
    # ---
    @property
    def pie(self):
        """
        The 'pie' property is a tuple of instances of
        Pie that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Pie
          - A list or tuple of dicts of string/value properties that
            will be passed to the Pie constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Pie]
        """
        return self["pie"]

    @pie.setter
    def pie(self, val):
        self["pie"] = val

    # pointcloud
    # ----------
    @property
    def pointcloud(self):
        """
        The 'pointcloud' property is a tuple of instances of
        Pointcloud that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Pointcloud
          - A list or tuple of dicts of string/value properties that
            will be passed to the Pointcloud constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Pointcloud]
        """
        return self["pointcloud"]

    @pointcloud.setter
    def pointcloud(self, val):
        self["pointcloud"] = val

    # sankey
    # ------
    @property
    def sankey(self):
        """
        The 'sankey' property is a tuple of instances of
        Sankey that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Sankey
          - A list or tuple of dicts of string/value properties that
            will be passed to the Sankey constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Sankey]
        """
        return self["sankey"]

    @sankey.setter
    def sankey(self, val):
        self["sankey"] = val

    # scatter3d
    # ---------
    @property
    def scatter3d(self):
        """
        The 'scatter3d' property is a tuple of instances of
        Scatter3d that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Scatter3d
          - A list or tuple of dicts of string/value properties that
            will be passed to the Scatter3d constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Scatter3d]
        """
        return self["scatter3d"]

    @scatter3d.setter
    def scatter3d(self, val):
        self["scatter3d"] = val

    # scattercarpet
    # -------------
    @property
    def scattercarpet(self):
        """
        The 'scattercarpet' property is a tuple of instances of
        Scattercarpet that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Scattercarpet
          - A list or tuple of dicts of string/value properties that
            will be passed to the Scattercarpet constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Scattercarpet]
        """
        return self["scattercarpet"]

    @scattercarpet.setter
    def scattercarpet(self, val):
        self["scattercarpet"] = val

    # scattergeo
    # ----------
    @property
    def scattergeo(self):
        """
        The 'scattergeo' property is a tuple of instances of
        Scattergeo that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Scattergeo
          - A list or tuple of dicts of string/value properties that
            will be passed to the Scattergeo constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Scattergeo]
        """
        return self["scattergeo"]

    @scattergeo.setter
    def scattergeo(self, val):
        self["scattergeo"] = val

    # scattergl
    # ---------
    @property
    def scattergl(self):
        """
        The 'scattergl' property is a tuple of instances of
        Scattergl that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Scattergl
          - A list or tuple of dicts of string/value properties that
            will be passed to the Scattergl constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Scattergl]
        """
        return self["scattergl"]

    @scattergl.setter
    def scattergl(self, val):
        self["scattergl"] = val

    # scattermapbox
    # -------------
    @property
    def scattermapbox(self):
        """
        The 'scattermapbox' property is a tuple of instances of
        Scattermapbox that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Scattermapbox
          - A list or tuple of dicts of string/value properties that
            will be passed to the Scattermapbox constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Scattermapbox]
        """
        return self["scattermapbox"]

    @scattermapbox.setter
    def scattermapbox(self, val):
        self["scattermapbox"] = val

    # scattermap
    # ----------
    @property
    def scattermap(self):
        """
        The 'scattermap' property is a tuple of instances of
        Scattermap that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Scattermap
          - A list or tuple of dicts of string/value properties that
            will be passed to the Scattermap constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Scattermap]
        """
        return self["scattermap"]

    @scattermap.setter
    def scattermap(self, val):
        self["scattermap"] = val

    # scatterpolargl
    # --------------
    @property
    def scatterpolargl(self):
        """
        The 'scatterpolargl' property is a tuple of instances of
        Scatterpolargl that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Scatterpolargl
          - A list or tuple of dicts of string/value properties that
            will be passed to the Scatterpolargl constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Scatterpolargl]
        """
        return self["scatterpolargl"]

    @scatterpolargl.setter
    def scatterpolargl(self, val):
        self["scatterpolargl"] = val

    # scatterpolar
    # ------------
    @property
    def scatterpolar(self):
        """
        The 'scatterpolar' property is a tuple of instances of
        Scatterpolar that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Scatterpolar
          - A list or tuple of dicts of string/value properties that
            will be passed to the Scatterpolar constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Scatterpolar]
        """
        return self["scatterpolar"]

    @scatterpolar.setter
    def scatterpolar(self, val):
        self["scatterpolar"] = val

    # scatter
    # -------
    @property
    def scatter(self):
        """
        The 'scatter' property is a tuple of instances of
        Scatter that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Scatter
          - A list or tuple of dicts of string/value properties that
            will be passed to the Scatter constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Scatter]
        """
        return self["scatter"]

    @scatter.setter
    def scatter(self, val):
        self["scatter"] = val

    # scattersmith
    # ------------
    @property
    def scattersmith(self):
        """
        The 'scattersmith' property is a tuple of instances of
        Scattersmith that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Scattersmith
          - A list or tuple of dicts of string/value properties that
            will be passed to the Scattersmith constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Scattersmith]
        """
        return self["scattersmith"]

    @scattersmith.setter
    def scattersmith(self, val):
        self["scattersmith"] = val

    # scatterternary
    # --------------
    @property
    def scatterternary(self):
        """
        The 'scatterternary' property is a tuple of instances of
        Scatterternary that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Scatterternary
          - A list or tuple of dicts of string/value properties that
            will be passed to the Scatterternary constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Scatterternary]
        """
        return self["scatterternary"]

    @scatterternary.setter
    def scatterternary(self, val):
        self["scatterternary"] = val

    # splom
    # -----
    @property
    def splom(self):
        """
        The 'splom' property is a tuple of instances of
        Splom that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Splom
          - A list or tuple of dicts of string/value properties that
            will be passed to the Splom constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Splom]
        """
        return self["splom"]

    @splom.setter
    def splom(self, val):
        self["splom"] = val

    # streamtube
    # ----------
    @property
    def streamtube(self):
        """
        The 'streamtube' property is a tuple of instances of
        Streamtube that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Streamtube
          - A list or tuple of dicts of string/value properties that
            will be passed to the Streamtube constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Streamtube]
        """
        return self["streamtube"]

    @streamtube.setter
    def streamtube(self, val):
        self["streamtube"] = val

    # sunburst
    # --------
    @property
    def sunburst(self):
        """
        The 'sunburst' property is a tuple of instances of
        Sunburst that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Sunburst
          - A list or tuple of dicts of string/value properties that
            will be passed to the Sunburst constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Sunburst]
        """
        return self["sunburst"]

    @sunburst.setter
    def sunburst(self, val):
        self["sunburst"] = val

    # surface
    # -------
    @property
    def surface(self):
        """
        The 'surface' property is a tuple of instances of
        Surface that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Surface
          - A list or tuple of dicts of string/value properties that
            will be passed to the Surface constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Surface]
        """
        return self["surface"]

    @surface.setter
    def surface(self, val):
        self["surface"] = val

    # table
    # -----
    @property
    def table(self):
        """
        The 'table' property is a tuple of instances of
        Table that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Table
          - A list or tuple of dicts of string/value properties that
            will be passed to the Table constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Table]
        """
        return self["table"]

    @table.setter
    def table(self, val):
        self["table"] = val

    # treemap
    # -------
    @property
    def treemap(self):
        """
        The 'treemap' property is a tuple of instances of
        Treemap that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Treemap
          - A list or tuple of dicts of string/value properties that
            will be passed to the Treemap constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Treemap]
        """
        return self["treemap"]

    @treemap.setter
    def treemap(self, val):
        self["treemap"] = val

    # violin
    # ------
    @property
    def violin(self):
        """
        The 'violin' property is a tuple of instances of
        Violin that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Violin
          - A list or tuple of dicts of string/value properties that
            will be passed to the Violin constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Violin]
        """
        return self["violin"]

    @violin.setter
    def violin(self, val):
        self["violin"] = val

    # volume
    # ------
    @property
    def volume(self):
        """
        The 'volume' property is a tuple of instances of
        Volume that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Volume
          - A list or tuple of dicts of string/value properties that
            will be passed to the Volume constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Volume]
        """
        return self["volume"]

    @volume.setter
    def volume(self, val):
        self["volume"] = val

    # waterfall
    # ---------
    @property
    def waterfall(self):
        """
        The 'waterfall' property is a tuple of instances of
        Waterfall that may be specified as:
          - A list or tuple of instances of plotly.graph_objs.layout.template.data.Waterfall
          - A list or tuple of dicts of string/value properties that
            will be passed to the Waterfall constructor

            Supported dict properties:

        Returns
        -------
        tuple[plotly.graph_objs.layout.template.data.Waterfall]
        """
        return self["waterfall"]

    @waterfall.setter
    def waterfall(self, val):
        self["waterfall"] = val

    # Self properties description
    # ---------------------------
    @property
    def _prop_descriptions(self):
        return """\
        barpolar
            A tuple of :class:`plotly.graph_objects.Barpolar`
            instances or dicts with compatible properties
        bar
            A tuple of :class:`plotly.graph_objects.Bar` instances
            or dicts with compatible properties
        box
            A tuple of :class:`plotly.graph_objects.Box` instances
            or dicts with compatible properties
        candlestick
            A tuple of :class:`plotly.graph_objects.Candlestick`
            instances or dicts with compatible properties
        carpet
            A tuple of :class:`plotly.graph_objects.Carpet`
            instances or dicts with compatible properties
        choroplethmapbox
            A tuple of
            :class:`plotly.graph_objects.Choroplethmapbox`
            instances or dicts with compatible properties
        choroplethmap
            A tuple of :class:`plotly.graph_objects.Choroplethmap`
            instances or dicts with compatible properties
        choropleth
            A tuple of :class:`plotly.graph_objects.Choropleth`
            instances or dicts with compatible properties
        cone
            A tuple of :class:`plotly.graph_objects.Cone` instances
            or dicts with compatible properties
        contourcarpet
            A tuple of :class:`plotly.graph_objects.Contourcarpet`
            instances or dicts with compatible properties
        contour
            A tuple of :class:`plotly.graph_objects.Contour`
            instances or dicts with compatible properties
        densitymapbox
            A tuple of :class:`plotly.graph_objects.Densitymapbox`
            instances or dicts with compatible properties
        densitymap
            A tuple of :class:`plotly.graph_objects.Densitymap`
            instances or dicts with compatible properties
        funnelarea
            A tuple of :class:`plotly.graph_objects.Funnelarea`
            instances or dicts with compatible properties
        funnel
            A tuple of :class:`plotly.graph_objects.Funnel`
            instances or dicts with compatible properties
        heatmapgl
            A tuple of :class:`plotly.graph_objects.Heatmapgl`
            instances or dicts with compatible properties
        heatmap
            A tuple of :class:`plotly.graph_objects.Heatmap`
            instances or dicts with compatible properties
        histogram2dcontour
            A tuple of
            :class:`plotly.graph_objects.Histogram2dContour`
            instances or dicts with compatible properties
        histogram2d
            A tuple of :class:`plotly.graph_objects.Histogram2d`
            instances or dicts with compatible properties
        histogram
            A tuple of :class:`plotly.graph_objects.Histogram`
            instances or dicts with compatible properties
        icicle
            A tuple of :class:`plotly.graph_objects.Icicle`
            instances or dicts with compatible properties
        image
            A tuple of :class:`plotly.graph_objects.Image`
            instances or dicts with compatible properties
        indicator
            A tuple of :class:`plotly.graph_objects.Indicator`
            instances or dicts with compatible properties
        isosurface
            A tuple of :class:`plotly.graph_objects.Isosurface`
            instances or dicts with compatible properties
        mesh3d
            A tuple of :class:`plotly.graph_objects.Mesh3d`
            instances or dicts with compatible properties
        ohlc
            A tuple of :class:`plotly.graph_objects.Ohlc` instances
            or dicts with compatible properties
        parcats
            A tuple of :class:`plotly.graph_objects.Parcats`
            instances or dicts with compatible properties
        parcoords
            A tuple of :class:`plotly.graph_objects.Parcoords`
            instances or dicts with compatible properties
        pie
            A tuple of :class:`plotly.graph_objects.Pie` instances
            or dicts with compatible properties
        pointcloud
            A tuple of :class:`plotly.graph_objects.Pointcloud`
            instances or dicts with compatible properties
        sankey
            A tuple of :class:`plotly.graph_objects.Sankey`
            instances or dicts with compatible properties
        scatter3d
            A tuple of :class:`plotly.graph_objects.Scatter3d`
            instances or dicts with compatible properties
        scattercarpet
            A tuple of :class:`plotly.graph_objects.Scattercarpet`
            instances or dicts with compatible properties
        scattergeo
            A tuple of :class:`plotly.graph_objects.Scattergeo`
            instances or dicts with compatible properties
        scattergl
            A tuple of :class:`plotly.graph_objects.Scattergl`
            instances or dicts with compatible properties
        scattermapbox
            A tuple of :class:`plotly.graph_objects.Scattermapbox`
            instances or dicts with compatible properties
        scattermap
            A tuple of :class:`plotly.graph_objects.Scattermap`
            instances or dicts with compatible properties
        scatterpolargl
            A tuple of :class:`plotly.graph_objects.Scatterpolargl`
            instances or dicts with compatible properties
        scatterpolar
            A tuple of :class:`plotly.graph_objects.Scatterpolar`
            instances or dicts with compatible properties
        scatter
            A tuple of :class:`plotly.graph_objects.Scatter`
            instances or dicts with compatible properties
        scattersmith
            A tuple of :class:`plotly.graph_objects.Scattersmith`
            instances or dicts with compatible properties
        scatterternary
            A tuple of :class:`plotly.graph_objects.Scatterternary`
            instances or dicts with compatible properties
        splom
            A tuple of :class:`plotly.graph_objects.Splom`
            instances or dicts with compatible properties
        streamtube
            A tuple of :class:`plotly.graph_objects.Streamtube`
            instances or dicts with compatible properties
        sunburst
            A tuple of :class:`plotly.graph_objects.Sunburst`
            instances or dicts with compatible properties
        surface
            A tuple of :class:`plotly.graph_objects.Surface`
            instances or dicts with compatible properties
        table
            A tuple of :class:`plotly.graph_objects.Table`
            instances or dicts with compatible properties
        treemap
            A tuple of :class:`plotly.graph_objects.Treemap`
            instances or dicts with compatible properties
        violin
            A tuple of :class:`plotly.graph_objects.Violin`
            instances or dicts with compatible properties
        volume
            A tuple of :class:`plotly.graph_objects.Volume`
            instances or dicts with compatible properties
        waterfall
            A tuple of :class:`plotly.graph_objects.Waterfall`
            instances or dicts with compatible properties
        """

    def __init__(
        self,
        arg=None,
        barpolar=None,
        bar=None,
        box=None,
        candlestick=None,
        carpet=None,
        choroplethmapbox=None,
        choroplethmap=None,
        choropleth=None,
        cone=None,
        contourcarpet=None,
        contour=None,
        densitymapbox=None,
        densitymap=None,
        funnelarea=None,
        funnel=None,
        heatmapgl=None,
        heatmap=None,
        histogram2dcontour=None,
        histogram2d=None,
        histogram=None,
        icicle=None,
        image=None,
        indicator=None,
        isosurface=None,
        mesh3d=None,
        ohlc=None,
        parcats=None,
        parcoords=None,
        pie=None,
        pointcloud=None,
        sankey=None,
        scatter3d=None,
        scattercarpet=None,
        scattergeo=None,
        scattergl=None,
        scattermapbox=None,
        scattermap=None,
        scatterpolargl=None,
        scatterpolar=None,
        scatter=None,
        scattersmith=None,
        scatterternary=None,
        splom=None,
        streamtube=None,
        sunburst=None,
        surface=None,
        table=None,
        treemap=None,
        violin=None,
        volume=None,
        waterfall=None,
        **kwargs,
    ):
        """
        Construct a new Data object

        Parameters
        ----------
        arg
            dict of properties compatible with this constructor or
            an instance of
            :class:`plotly.graph_objs.layout.template.Data`
        barpolar
            A tuple of :class:`plotly.graph_objects.Barpolar`
            instances or dicts with compatible properties
        bar
            A tuple of :class:`plotly.graph_objects.Bar` instances
            or dicts with compatible properties
        box
            A tuple of :class:`plotly.graph_objects.Box` instances
            or dicts with compatible properties
        candlestick
            A tuple of :class:`plotly.graph_objects.Candlestick`
            instances or dicts with compatible properties
        carpet
            A tuple of :class:`plotly.graph_objects.Carpet`
            instances or dicts with compatible properties
        choroplethmapbox
            A tuple of
            :class:`plotly.graph_objects.Choroplethmapbox`
            instances or dicts with compatible properties
        choroplethmap
            A tuple of :class:`plotly.graph_objects.Choroplethmap`
            instances or dicts with compatible properties
        choropleth
            A tuple of :class:`plotly.graph_objects.Choropleth`
            instances or dicts with compatible properties
        cone
            A tuple of :class:`plotly.graph_objects.Cone` instances
            or dicts with compatible properties
        contourcarpet
            A tuple of :class:`plotly.graph_objects.Contourcarpet`
            instances or dicts with compatible properties
        contour
            A tuple of :class:`plotly.graph_objects.Contour`
            instances or dicts with compatible properties
        densitymapbox
            A tuple of :class:`plotly.graph_objects.Densitymapbox`
            instances or dicts with compatible properties
        densitymap
            A tuple of :class:`plotly.graph_objects.Densitymap`
            instances or dicts with compatible properties
        funnelarea
            A tuple of :class:`plotly.graph_objects.Funnelarea`
            instances or dicts with compatible properties
        funnel
            A tuple of :class:`plotly.graph_objects.Funnel`
            instances or dicts with compatible properties
        heatmapgl
            A tuple of :class:`plotly.graph_objects.Heatmapgl`
            instances or dicts with compatible properties
        heatmap
            A tuple of :class:`plotly.graph_objects.Heatmap`
            instances or dicts with compatible properties
        histogram2dcontour
            A tuple of
            :class:`plotly.graph_objects.Histogram2dContour`
            instances or dicts with compatible properties
        histogram2d
            A tuple of :class:`plotly.graph_objects.Histogram2d`
            instances or dicts with compatible properties
        histogram
            A tuple of :class:`plotly.graph_objects.Histogram`
            instances or dicts with compatible properties
        icicle
            A tuple of :class:`plotly.graph_objects.Icicle`
            instances or dicts with compatible properties
        image
            A tuple of :class:`plotly.graph_objects.Image`
            instances or dicts with compatible properties
        indicator
            A tuple of :class:`plotly.graph_objects.Indicator`
            instances or dicts with compatible properties
        isosurface
            A tuple of :class:`plotly.graph_objects.Isosurface`
            instances or dicts with compatible properties
        mesh3d
            A tuple of :class:`plotly.graph_objects.Mesh3d`
            instances or dicts with compatible properties
        ohlc
            A tuple of :class:`plotly.graph_objects.Ohlc` instances
            or dicts with compatible properties
        parcats
            A tuple of :class:`plotly.graph_objects.Parcats`
            instances or dicts with compatible properties
        parcoords
            A tuple of :class:`plotly.graph_objects.Parcoords`
            instances or dicts with compatible properties
        pie
            A tuple of :class:`plotly.graph_objects.Pie` instances
            or dicts with compatible properties
        pointcloud
            A tuple of :class:`plotly.graph_objects.Pointcloud`
            instances or dicts with compatible properties
        sankey
            A tuple of :class:`plotly.graph_objects.Sankey`
            instances or dicts with compatible properties
        scatter3d
            A tuple of :class:`plotly.graph_objects.Scatter3d`
            instances or dicts with compatible properties
        scattercarpet
            A tuple of :class:`plotly.graph_objects.Scattercarpet`
            instances or dicts with compatible properties
        scattergeo
            A tuple of :class:`plotly.graph_objects.Scattergeo`
            instances or dicts with compatible properties
        scattergl
            A tuple of :class:`plotly.graph_objects.Scattergl`
            instances or dicts with compatible properties
        scattermapbox
            A tuple of :class:`plotly.graph_objects.Scattermapbox`
            instances or dicts with compatible properties
        scattermap
            A tuple of :class:`plotly.graph_objects.Scattermap`
            instances or dicts with compatible properties
        scatterpolargl
            A tuple of :class:`plotly.graph_objects.Scatterpolargl`
            instances or dicts with compatible properties
        scatterpolar
            A tuple of :class:`plotly.graph_objects.Scatterpolar`
            instances or dicts with compatible properties
        scatter
            A tuple of :class:`plotly.graph_objects.Scatter`
            instances or dicts with compatible properties
        scattersmith
            A tuple of :class:`plotly.graph_objects.Scattersmith`
            instances or dicts with compatible properties
        scatterternary
            A tuple of :class:`plotly.graph_objects.Scatterternary`
            instances or dicts with compatible properties
        splom
            A tuple of :class:`plotly.graph_objects.Splom`
            instances or dicts with compatible properties
        streamtube
            A tuple of :class:`plotly.graph_objects.Streamtube`
            instances or dicts with compatible properties
        sunburst
            A tuple of :class:`plotly.graph_objects.Sunburst`
            instances or dicts with compatible properties
        surface
            A tuple of :class:`plotly.graph_objects.Surface`
            instances or dicts with compatible properties
        table
            A tuple of :class:`plotly.graph_objects.Table`
            instances or dicts with compatible properties
        treemap
            A tuple of :class:`plotly.graph_objects.Treemap`
            instances or dicts with compatible properties
        violin
            A tuple of :class:`plotly.graph_objects.Violin`
            instances or dicts with compatible properties
        volume
            A tuple of :class:`plotly.graph_objects.Volume`
            instances or dicts with compatible properties
        waterfall
            A tuple of :class:`plotly.graph_objects.Waterfall`
            instances or dicts with compatible properties

        Returns
        -------
        Data
        """
        super(Data, self).__init__("data")

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
The first argument to the plotly.graph_objs.layout.template.Data
constructor must be a dict or
an instance of :class:`plotly.graph_objs.layout.template.Data`"""
            )

        # Handle skip_invalid
        # -------------------
        self._skip_invalid = kwargs.pop("skip_invalid", False)
        self._validate = kwargs.pop("_validate", True)

        # Populate data dict with properties
        # ----------------------------------
        _v = arg.pop("barpolar", None)
        _v = barpolar if barpolar is not None else _v
        if _v is not None:
            self["barpolar"] = _v
        _v = arg.pop("bar", None)
        _v = bar if bar is not None else _v
        if _v is not None:
            self["bar"] = _v
        _v = arg.pop("box", None)
        _v = box if box is not None else _v
        if _v is not None:
            self["box"] = _v
        _v = arg.pop("candlestick", None)
        _v = candlestick if candlestick is not None else _v
        if _v is not None:
            self["candlestick"] = _v
        _v = arg.pop("carpet", None)
        _v = carpet if carpet is not None else _v
        if _v is not None:
            self["carpet"] = _v
        _v = arg.pop("choroplethmapbox", None)
        _v = choroplethmapbox if choroplethmapbox is not None else _v
        if _v is not None:
            self["choroplethmapbox"] = _v
        _v = arg.pop("choroplethmap", None)
        _v = choroplethmap if choroplethmap is not None else _v
        if _v is not None:
            self["choroplethmap"] = _v
        _v = arg.pop("choropleth", None)
        _v = choropleth if choropleth is not None else _v
        if _v is not None:
            self["choropleth"] = _v
        _v = arg.pop("cone", None)
        _v = cone if cone is not None else _v
        if _v is not None:
            self["cone"] = _v
        _v = arg.pop("contourcarpet", None)
        _v = contourcarpet if contourcarpet is not None else _v
        if _v is not None:
            self["contourcarpet"] = _v
        _v = arg.pop("contour", None)
        _v = contour if contour is not None else _v
        if _v is not None:
            self["contour"] = _v
        _v = arg.pop("densitymapbox", None)
        _v = densitymapbox if densitymapbox is not None else _v
        if _v is not None:
            self["densitymapbox"] = _v
        _v = arg.pop("densitymap", None)
        _v = densitymap if densitymap is not None else _v
        if _v is not None:
            self["densitymap"] = _v
        _v = arg.pop("funnelarea", None)
        _v = funnelarea if funnelarea is not None else _v
        if _v is not None:
            self["funnelarea"] = _v
        _v = arg.pop("funnel", None)
        _v = funnel if funnel is not None else _v
        if _v is not None:
            self["funnel"] = _v
        _v = arg.pop("heatmapgl", None)
        _v = heatmapgl if heatmapgl is not None else _v
        if _v is not None:
            self["heatmapgl"] = _v
        _v = arg.pop("heatmap", None)
        _v = heatmap if heatmap is not None else _v
        if _v is not None:
            self["heatmap"] = _v
        _v = arg.pop("histogram2dcontour", None)
        _v = histogram2dcontour if histogram2dcontour is not None else _v
        if _v is not None:
            self["histogram2dcontour"] = _v
        _v = arg.pop("histogram2d", None)
        _v = histogram2d if histogram2d is not None else _v
        if _v is not None:
            self["histogram2d"] = _v
        _v = arg.pop("histogram", None)
        _v = histogram if histogram is not None else _v
        if _v is not None:
            self["histogram"] = _v
        _v = arg.pop("icicle", None)
        _v = icicle if icicle is not None else _v
        if _v is not None:
            self["icicle"] = _v
        _v = arg.pop("image", None)
        _v = image if image is not None else _v
        if _v is not None:
            self["image"] = _v
        _v = arg.pop("indicator", None)
        _v = indicator if indicator is not None else _v
        if _v is not None:
            self["indicator"] = _v
        _v = arg.pop("isosurface", None)
        _v = isosurface if isosurface is not None else _v
        if _v is not None:
            self["isosurface"] = _v
        _v = arg.pop("mesh3d", None)
        _v = mesh3d if mesh3d is not None else _v
        if _v is not None:
            self["mesh3d"] = _v
        _v = arg.pop("ohlc", None)
        _v = ohlc if ohlc is not None else _v
        if _v is not None:
            self["ohlc"] = _v
        _v = arg.pop("parcats", None)
        _v = parcats if parcats is not None else _v
        if _v is not None:
            self["parcats"] = _v
        _v = arg.pop("parcoords", None)
        _v = parcoords if parcoords is not None else _v
        if _v is not None:
            self["parcoords"] = _v
        _v = arg.pop("pie", None)
        _v = pie if pie is not None else _v
        if _v is not None:
            self["pie"] = _v
        _v = arg.pop("pointcloud", None)
        _v = pointcloud if pointcloud is not None else _v
        if _v is not None:
            self["pointcloud"] = _v
        _v = arg.pop("sankey", None)
        _v = sankey if sankey is not None else _v
        if _v is not None:
            self["sankey"] = _v
        _v = arg.pop("scatter3d", None)
        _v = scatter3d if scatter3d is not None else _v
        if _v is not None:
            self["scatter3d"] = _v
        _v = arg.pop("scattercarpet", None)
        _v = scattercarpet if scattercarpet is not None else _v
        if _v is not None:
            self["scattercarpet"] = _v
        _v = arg.pop("scattergeo", None)
        _v = scattergeo if scattergeo is not None else _v
        if _v is not None:
            self["scattergeo"] = _v
        _v = arg.pop("scattergl", None)
        _v = scattergl if scattergl is not None else _v
        if _v is not None:
            self["scattergl"] = _v
        _v = arg.pop("scattermapbox", None)
        _v = scattermapbox if scattermapbox is not None else _v
        if _v is not None:
            self["scattermapbox"] = _v
        _v = arg.pop("scattermap", None)
        _v = scattermap if scattermap is not None else _v
        if _v is not None:
            self["scattermap"] = _v
        _v = arg.pop("scatterpolargl", None)
        _v = scatterpolargl if scatterpolargl is not None else _v
        if _v is not None:
            self["scatterpolargl"] = _v
        _v = arg.pop("scatterpolar", None)
        _v = scatterpolar if scatterpolar is not None else _v
        if _v is not None:
            self["scatterpolar"] = _v
        _v = arg.pop("scatter", None)
        _v = scatter if scatter is not None else _v
        if _v is not None:
            self["scatter"] = _v
        _v = arg.pop("scattersmith", None)
        _v = scattersmith if scattersmith is not None else _v
        if _v is not None:
            self["scattersmith"] = _v
        _v = arg.pop("scatterternary", None)
        _v = scatterternary if scatterternary is not None else _v
        if _v is not None:
            self["scatterternary"] = _v
        _v = arg.pop("splom", None)
        _v = splom if splom is not None else _v
        if _v is not None:
            self["splom"] = _v
        _v = arg.pop("streamtube", None)
        _v = streamtube if streamtube is not None else _v
        if _v is not None:
            self["streamtube"] = _v
        _v = arg.pop("sunburst", None)
        _v = sunburst if sunburst is not None else _v
        if _v is not None:
            self["sunburst"] = _v
        _v = arg.pop("surface", None)
        _v = surface if surface is not None else _v
        if _v is not None:
            self["surface"] = _v
        _v = arg.pop("table", None)
        _v = table if table is not None else _v
        if _v is not None:
            self["table"] = _v
        _v = arg.pop("treemap", None)
        _v = treemap if treemap is not None else _v
        if _v is not None:
            self["treemap"] = _v
        _v = arg.pop("violin", None)
        _v = violin if violin is not None else _v
        if _v is not None:
            self["violin"] = _v
        _v = arg.pop("volume", None)
        _v = volume if volume is not None else _v
        if _v is not None:
            self["volume"] = _v
        _v = arg.pop("waterfall", None)
        _v = waterfall if waterfall is not None else _v
        if _v is not None:
            self["waterfall"] = _v

        # Process unknown kwargs
        # ----------------------
        self._process_kwargs(**dict(arg, **kwargs))

        # Reset skip_invalid
        # ------------------
        self._skip_invalid = False
