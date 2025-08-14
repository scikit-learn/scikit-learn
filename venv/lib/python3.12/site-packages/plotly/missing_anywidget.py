from .basedatatypes import BaseFigure


class FigureWidget(BaseFigure):
    """
    FigureWidget stand-in for use when anywidget is not installed. The only purpose
    of this class is to provide something to import as
    `plotly.graph_objs.FigureWidget` when anywidget is not installed. This class
    simply raises an informative error message when the constructor is called
    """

    def __init__(self, *args, **kwargs):
        raise ImportError("Please install anywidget to use the FigureWidget class")
