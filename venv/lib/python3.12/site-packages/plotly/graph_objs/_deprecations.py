import warnings

warnings.filterwarnings(
    "default", r"plotly\.graph_objs\.\w+ is deprecated", DeprecationWarning
)


class Data(list):
    """
        plotly.graph_objs.Data is deprecated.
    Please replace it with a list or tuple of instances of the following types
      - plotly.graph_objs.Scatter
      - plotly.graph_objs.Bar
      - plotly.graph_objs.Area
      - plotly.graph_objs.Histogram
      - etc.

    """

    def __init__(self, *args, **kwargs):
        """
                plotly.graph_objs.Data is deprecated.
        Please replace it with a list or tuple of instances of the following types
          - plotly.graph_objs.Scatter
          - plotly.graph_objs.Bar
          - plotly.graph_objs.Area
          - plotly.graph_objs.Histogram
          - etc.

        """
        warnings.warn(
            """plotly.graph_objs.Data is deprecated.
Please replace it with a list or tuple of instances of the following types
  - plotly.graph_objs.Scatter
  - plotly.graph_objs.Bar
  - plotly.graph_objs.Area
  - plotly.graph_objs.Histogram
  - etc.
""",
            DeprecationWarning,
        )
        super().__init__(*args, **kwargs)


class Annotations(list):
    """
        plotly.graph_objs.Annotations is deprecated.
    Please replace it with a list or tuple of instances of the following types
      - plotly.graph_objs.layout.Annotation
      - plotly.graph_objs.layout.scene.Annotation

    """

    def __init__(self, *args, **kwargs):
        """
                plotly.graph_objs.Annotations is deprecated.
        Please replace it with a list or tuple of instances of the following types
          - plotly.graph_objs.layout.Annotation
          - plotly.graph_objs.layout.scene.Annotation

        """
        warnings.warn(
            """plotly.graph_objs.Annotations is deprecated.
Please replace it with a list or tuple of instances of the following types
  - plotly.graph_objs.layout.Annotation
  - plotly.graph_objs.layout.scene.Annotation
""",
            DeprecationWarning,
        )
        super().__init__(*args, **kwargs)


class Frames(list):
    """
        plotly.graph_objs.Frames is deprecated.
    Please replace it with a list or tuple of instances of the following types
      - plotly.graph_objs.Frame

    """

    def __init__(self, *args, **kwargs):
        """
                plotly.graph_objs.Frames is deprecated.
        Please replace it with a list or tuple of instances of the following types
          - plotly.graph_objs.Frame

        """
        warnings.warn(
            """plotly.graph_objs.Frames is deprecated.
Please replace it with a list or tuple of instances of the following types
  - plotly.graph_objs.Frame
""",
            DeprecationWarning,
        )
        super().__init__(*args, **kwargs)


class AngularAxis(dict):
    """
        plotly.graph_objs.AngularAxis is deprecated.
    Please replace it with one of the following more specific types
      - plotly.graph_objs.layout.AngularAxis
      - plotly.graph_objs.layout.polar.AngularAxis

    """

    def __init__(self, *args, **kwargs):
        """
                plotly.graph_objs.AngularAxis is deprecated.
        Please replace it with one of the following more specific types
          - plotly.graph_objs.layout.AngularAxis
          - plotly.graph_objs.layout.polar.AngularAxis

        """
        warnings.warn(
            """plotly.graph_objs.AngularAxis is deprecated.
Please replace it with one of the following more specific types
  - plotly.graph_objs.layout.AngularAxis
  - plotly.graph_objs.layout.polar.AngularAxis
""",
            DeprecationWarning,
        )
        super().__init__(*args, **kwargs)


class Annotation(dict):
    """
        plotly.graph_objs.Annotation is deprecated.
    Please replace it with one of the following more specific types
      - plotly.graph_objs.layout.Annotation
      - plotly.graph_objs.layout.scene.Annotation

    """

    def __init__(self, *args, **kwargs):
        """
                plotly.graph_objs.Annotation is deprecated.
        Please replace it with one of the following more specific types
          - plotly.graph_objs.layout.Annotation
          - plotly.graph_objs.layout.scene.Annotation

        """
        warnings.warn(
            """plotly.graph_objs.Annotation is deprecated.
Please replace it with one of the following more specific types
  - plotly.graph_objs.layout.Annotation
  - plotly.graph_objs.layout.scene.Annotation
""",
            DeprecationWarning,
        )
        super().__init__(*args, **kwargs)


class ColorBar(dict):
    """
        plotly.graph_objs.ColorBar is deprecated.
    Please replace it with one of the following more specific types
      - plotly.graph_objs.scatter.marker.ColorBar
      - plotly.graph_objs.surface.ColorBar
      - etc.

    """

    def __init__(self, *args, **kwargs):
        """
                plotly.graph_objs.ColorBar is deprecated.
        Please replace it with one of the following more specific types
          - plotly.graph_objs.scatter.marker.ColorBar
          - plotly.graph_objs.surface.ColorBar
          - etc.

        """
        warnings.warn(
            """plotly.graph_objs.ColorBar is deprecated.
Please replace it with one of the following more specific types
  - plotly.graph_objs.scatter.marker.ColorBar
  - plotly.graph_objs.surface.ColorBar
  - etc.
""",
            DeprecationWarning,
        )
        super().__init__(*args, **kwargs)


class Contours(dict):
    """
        plotly.graph_objs.Contours is deprecated.
    Please replace it with one of the following more specific types
      - plotly.graph_objs.contour.Contours
      - plotly.graph_objs.surface.Contours
      - etc.

    """

    def __init__(self, *args, **kwargs):
        """
                plotly.graph_objs.Contours is deprecated.
        Please replace it with one of the following more specific types
          - plotly.graph_objs.contour.Contours
          - plotly.graph_objs.surface.Contours
          - etc.

        """
        warnings.warn(
            """plotly.graph_objs.Contours is deprecated.
Please replace it with one of the following more specific types
  - plotly.graph_objs.contour.Contours
  - plotly.graph_objs.surface.Contours
  - etc.
""",
            DeprecationWarning,
        )
        super().__init__(*args, **kwargs)


class ErrorX(dict):
    """
        plotly.graph_objs.ErrorX is deprecated.
    Please replace it with one of the following more specific types
      - plotly.graph_objs.scatter.ErrorX
      - plotly.graph_objs.histogram.ErrorX
      - etc.

    """

    def __init__(self, *args, **kwargs):
        """
                plotly.graph_objs.ErrorX is deprecated.
        Please replace it with one of the following more specific types
          - plotly.graph_objs.scatter.ErrorX
          - plotly.graph_objs.histogram.ErrorX
          - etc.

        """
        warnings.warn(
            """plotly.graph_objs.ErrorX is deprecated.
Please replace it with one of the following more specific types
  - plotly.graph_objs.scatter.ErrorX
  - plotly.graph_objs.histogram.ErrorX
  - etc.
""",
            DeprecationWarning,
        )
        super().__init__(*args, **kwargs)


class ErrorY(dict):
    """
        plotly.graph_objs.ErrorY is deprecated.
    Please replace it with one of the following more specific types
      - plotly.graph_objs.scatter.ErrorY
      - plotly.graph_objs.histogram.ErrorY
      - etc.

    """

    def __init__(self, *args, **kwargs):
        """
                plotly.graph_objs.ErrorY is deprecated.
        Please replace it with one of the following more specific types
          - plotly.graph_objs.scatter.ErrorY
          - plotly.graph_objs.histogram.ErrorY
          - etc.

        """
        warnings.warn(
            """plotly.graph_objs.ErrorY is deprecated.
Please replace it with one of the following more specific types
  - plotly.graph_objs.scatter.ErrorY
  - plotly.graph_objs.histogram.ErrorY
  - etc.
""",
            DeprecationWarning,
        )
        super().__init__(*args, **kwargs)


class ErrorZ(dict):
    """
        plotly.graph_objs.ErrorZ is deprecated.
    Please replace it with one of the following more specific types
      - plotly.graph_objs.scatter3d.ErrorZ

    """

    def __init__(self, *args, **kwargs):
        """
                plotly.graph_objs.ErrorZ is deprecated.
        Please replace it with one of the following more specific types
          - plotly.graph_objs.scatter3d.ErrorZ

        """
        warnings.warn(
            """plotly.graph_objs.ErrorZ is deprecated.
Please replace it with one of the following more specific types
  - plotly.graph_objs.scatter3d.ErrorZ
""",
            DeprecationWarning,
        )
        super().__init__(*args, **kwargs)


class Font(dict):
    """
        plotly.graph_objs.Font is deprecated.
    Please replace it with one of the following more specific types
      - plotly.graph_objs.layout.Font
      - plotly.graph_objs.layout.hoverlabel.Font
      - etc.

    """

    def __init__(self, *args, **kwargs):
        """
                plotly.graph_objs.Font is deprecated.
        Please replace it with one of the following more specific types
          - plotly.graph_objs.layout.Font
          - plotly.graph_objs.layout.hoverlabel.Font
          - etc.

        """
        warnings.warn(
            """plotly.graph_objs.Font is deprecated.
Please replace it with one of the following more specific types
  - plotly.graph_objs.layout.Font
  - plotly.graph_objs.layout.hoverlabel.Font
  - etc.
""",
            DeprecationWarning,
        )
        super().__init__(*args, **kwargs)


class Legend(dict):
    """
        plotly.graph_objs.Legend is deprecated.
    Please replace it with one of the following more specific types
      - plotly.graph_objs.layout.Legend

    """

    def __init__(self, *args, **kwargs):
        """
                plotly.graph_objs.Legend is deprecated.
        Please replace it with one of the following more specific types
          - plotly.graph_objs.layout.Legend

        """
        warnings.warn(
            """plotly.graph_objs.Legend is deprecated.
Please replace it with one of the following more specific types
  - plotly.graph_objs.layout.Legend
""",
            DeprecationWarning,
        )
        super().__init__(*args, **kwargs)


class Line(dict):
    """
        plotly.graph_objs.Line is deprecated.
    Please replace it with one of the following more specific types
      - plotly.graph_objs.scatter.Line
      - plotly.graph_objs.layout.shape.Line
      - etc.

    """

    def __init__(self, *args, **kwargs):
        """
                plotly.graph_objs.Line is deprecated.
        Please replace it with one of the following more specific types
          - plotly.graph_objs.scatter.Line
          - plotly.graph_objs.layout.shape.Line
          - etc.

        """
        warnings.warn(
            """plotly.graph_objs.Line is deprecated.
Please replace it with one of the following more specific types
  - plotly.graph_objs.scatter.Line
  - plotly.graph_objs.layout.shape.Line
  - etc.
""",
            DeprecationWarning,
        )
        super().__init__(*args, **kwargs)


class Margin(dict):
    """
        plotly.graph_objs.Margin is deprecated.
    Please replace it with one of the following more specific types
      - plotly.graph_objs.layout.Margin

    """

    def __init__(self, *args, **kwargs):
        """
                plotly.graph_objs.Margin is deprecated.
        Please replace it with one of the following more specific types
          - plotly.graph_objs.layout.Margin

        """
        warnings.warn(
            """plotly.graph_objs.Margin is deprecated.
Please replace it with one of the following more specific types
  - plotly.graph_objs.layout.Margin
""",
            DeprecationWarning,
        )
        super().__init__(*args, **kwargs)


class Marker(dict):
    """
        plotly.graph_objs.Marker is deprecated.
    Please replace it with one of the following more specific types
      - plotly.graph_objs.scatter.Marker
      - plotly.graph_objs.histogram.selected.Marker
      - etc.

    """

    def __init__(self, *args, **kwargs):
        """
                plotly.graph_objs.Marker is deprecated.
        Please replace it with one of the following more specific types
          - plotly.graph_objs.scatter.Marker
          - plotly.graph_objs.histogram.selected.Marker
          - etc.

        """
        warnings.warn(
            """plotly.graph_objs.Marker is deprecated.
Please replace it with one of the following more specific types
  - plotly.graph_objs.scatter.Marker
  - plotly.graph_objs.histogram.selected.Marker
  - etc.
""",
            DeprecationWarning,
        )
        super().__init__(*args, **kwargs)


class RadialAxis(dict):
    """
        plotly.graph_objs.RadialAxis is deprecated.
    Please replace it with one of the following more specific types
      - plotly.graph_objs.layout.RadialAxis
      - plotly.graph_objs.layout.polar.RadialAxis

    """

    def __init__(self, *args, **kwargs):
        """
                plotly.graph_objs.RadialAxis is deprecated.
        Please replace it with one of the following more specific types
          - plotly.graph_objs.layout.RadialAxis
          - plotly.graph_objs.layout.polar.RadialAxis

        """
        warnings.warn(
            """plotly.graph_objs.RadialAxis is deprecated.
Please replace it with one of the following more specific types
  - plotly.graph_objs.layout.RadialAxis
  - plotly.graph_objs.layout.polar.RadialAxis
""",
            DeprecationWarning,
        )
        super().__init__(*args, **kwargs)


class Scene(dict):
    """
        plotly.graph_objs.Scene is deprecated.
    Please replace it with one of the following more specific types
      - plotly.graph_objs.layout.Scene

    """

    def __init__(self, *args, **kwargs):
        """
                plotly.graph_objs.Scene is deprecated.
        Please replace it with one of the following more specific types
          - plotly.graph_objs.layout.Scene

        """
        warnings.warn(
            """plotly.graph_objs.Scene is deprecated.
Please replace it with one of the following more specific types
  - plotly.graph_objs.layout.Scene
""",
            DeprecationWarning,
        )
        super().__init__(*args, **kwargs)


class Stream(dict):
    """
        plotly.graph_objs.Stream is deprecated.
    Please replace it with one of the following more specific types
      - plotly.graph_objs.scatter.Stream
      - plotly.graph_objs.area.Stream

    """

    def __init__(self, *args, **kwargs):
        """
                plotly.graph_objs.Stream is deprecated.
        Please replace it with one of the following more specific types
          - plotly.graph_objs.scatter.Stream
          - plotly.graph_objs.area.Stream

        """
        warnings.warn(
            """plotly.graph_objs.Stream is deprecated.
Please replace it with one of the following more specific types
  - plotly.graph_objs.scatter.Stream
  - plotly.graph_objs.area.Stream
""",
            DeprecationWarning,
        )
        super().__init__(*args, **kwargs)


class XAxis(dict):
    """
        plotly.graph_objs.XAxis is deprecated.
    Please replace it with one of the following more specific types
      - plotly.graph_objs.layout.XAxis
      - plotly.graph_objs.layout.scene.XAxis

    """

    def __init__(self, *args, **kwargs):
        """
                plotly.graph_objs.XAxis is deprecated.
        Please replace it with one of the following more specific types
          - plotly.graph_objs.layout.XAxis
          - plotly.graph_objs.layout.scene.XAxis

        """
        warnings.warn(
            """plotly.graph_objs.XAxis is deprecated.
Please replace it with one of the following more specific types
  - plotly.graph_objs.layout.XAxis
  - plotly.graph_objs.layout.scene.XAxis
""",
            DeprecationWarning,
        )
        super().__init__(*args, **kwargs)


class YAxis(dict):
    """
        plotly.graph_objs.YAxis is deprecated.
    Please replace it with one of the following more specific types
      - plotly.graph_objs.layout.YAxis
      - plotly.graph_objs.layout.scene.YAxis

    """

    def __init__(self, *args, **kwargs):
        """
                plotly.graph_objs.YAxis is deprecated.
        Please replace it with one of the following more specific types
          - plotly.graph_objs.layout.YAxis
          - plotly.graph_objs.layout.scene.YAxis

        """
        warnings.warn(
            """plotly.graph_objs.YAxis is deprecated.
Please replace it with one of the following more specific types
  - plotly.graph_objs.layout.YAxis
  - plotly.graph_objs.layout.scene.YAxis
""",
            DeprecationWarning,
        )
        super().__init__(*args, **kwargs)


class ZAxis(dict):
    """
        plotly.graph_objs.ZAxis is deprecated.
    Please replace it with one of the following more specific types
      - plotly.graph_objs.layout.scene.ZAxis

    """

    def __init__(self, *args, **kwargs):
        """
                plotly.graph_objs.ZAxis is deprecated.
        Please replace it with one of the following more specific types
          - plotly.graph_objs.layout.scene.ZAxis

        """
        warnings.warn(
            """plotly.graph_objs.ZAxis is deprecated.
Please replace it with one of the following more specific types
  - plotly.graph_objs.layout.scene.ZAxis
""",
            DeprecationWarning,
        )
        super().__init__(*args, **kwargs)


class XBins(dict):
    """
        plotly.graph_objs.XBins is deprecated.
    Please replace it with one of the following more specific types
      - plotly.graph_objs.histogram.XBins
      - plotly.graph_objs.histogram2d.XBins

    """

    def __init__(self, *args, **kwargs):
        """
                plotly.graph_objs.XBins is deprecated.
        Please replace it with one of the following more specific types
          - plotly.graph_objs.histogram.XBins
          - plotly.graph_objs.histogram2d.XBins

        """
        warnings.warn(
            """plotly.graph_objs.XBins is deprecated.
Please replace it with one of the following more specific types
  - plotly.graph_objs.histogram.XBins
  - plotly.graph_objs.histogram2d.XBins
""",
            DeprecationWarning,
        )
        super().__init__(*args, **kwargs)


class YBins(dict):
    """
        plotly.graph_objs.YBins is deprecated.
    Please replace it with one of the following more specific types
      - plotly.graph_objs.histogram.YBins
      - plotly.graph_objs.histogram2d.YBins

    """

    def __init__(self, *args, **kwargs):
        """
                plotly.graph_objs.YBins is deprecated.
        Please replace it with one of the following more specific types
          - plotly.graph_objs.histogram.YBins
          - plotly.graph_objs.histogram2d.YBins

        """
        warnings.warn(
            """plotly.graph_objs.YBins is deprecated.
Please replace it with one of the following more specific types
  - plotly.graph_objs.histogram.YBins
  - plotly.graph_objs.histogram2d.YBins
""",
            DeprecationWarning,
        )
        super().__init__(*args, **kwargs)


class Trace(dict):
    """
        plotly.graph_objs.Trace is deprecated.
    Please replace it with one of the following more specific types
      - plotly.graph_objs.Scatter
      - plotly.graph_objs.Bar
      - plotly.graph_objs.Area
      - plotly.graph_objs.Histogram
      - etc.

    """

    def __init__(self, *args, **kwargs):
        """
                plotly.graph_objs.Trace is deprecated.
        Please replace it with one of the following more specific types
          - plotly.graph_objs.Scatter
          - plotly.graph_objs.Bar
          - plotly.graph_objs.Area
          - plotly.graph_objs.Histogram
          - etc.

        """
        warnings.warn(
            """plotly.graph_objs.Trace is deprecated.
Please replace it with one of the following more specific types
  - plotly.graph_objs.Scatter
  - plotly.graph_objs.Bar
  - plotly.graph_objs.Area
  - plotly.graph_objs.Histogram
  - etc.
""",
            DeprecationWarning,
        )
        super().__init__(*args, **kwargs)


class Histogram2dcontour(dict):
    """
        plotly.graph_objs.Histogram2dcontour is deprecated.
    Please replace it with one of the following more specific types
      - plotly.graph_objs.Histogram2dContour

    """

    def __init__(self, *args, **kwargs):
        """
                plotly.graph_objs.Histogram2dcontour is deprecated.
        Please replace it with one of the following more specific types
          - plotly.graph_objs.Histogram2dContour

        """
        warnings.warn(
            """plotly.graph_objs.Histogram2dcontour is deprecated.
Please replace it with one of the following more specific types
  - plotly.graph_objs.Histogram2dContour
""",
            DeprecationWarning,
        )
        super().__init__(*args, **kwargs)
