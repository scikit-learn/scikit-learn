import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._activeselection import Activeselection
    from ._activeshape import Activeshape
    from ._annotation import Annotation
    from ._coloraxis import Coloraxis
    from ._colorscale import Colorscale
    from ._font import Font
    from ._geo import Geo
    from ._grid import Grid
    from ._hoverlabel import Hoverlabel
    from ._image import Image
    from ._legend import Legend
    from ._map import Map
    from ._mapbox import Mapbox
    from ._margin import Margin
    from ._modebar import Modebar
    from ._newselection import Newselection
    from ._newshape import Newshape
    from ._polar import Polar
    from ._scene import Scene
    from ._selection import Selection
    from ._shape import Shape
    from ._slider import Slider
    from ._smith import Smith
    from ._template import Template
    from ._ternary import Ternary
    from ._title import Title
    from ._transition import Transition
    from ._uniformtext import Uniformtext
    from ._updatemenu import Updatemenu
    from ._xaxis import XAxis
    from ._yaxis import YAxis
    from . import annotation
    from . import coloraxis
    from . import geo
    from . import grid
    from . import hoverlabel
    from . import legend
    from . import map
    from . import mapbox
    from . import newselection
    from . import newshape
    from . import polar
    from . import scene
    from . import selection
    from . import shape
    from . import slider
    from . import smith
    from . import template
    from . import ternary
    from . import title
    from . import updatemenu
    from . import xaxis
    from . import yaxis
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [
            ".annotation",
            ".coloraxis",
            ".geo",
            ".grid",
            ".hoverlabel",
            ".legend",
            ".map",
            ".mapbox",
            ".newselection",
            ".newshape",
            ".polar",
            ".scene",
            ".selection",
            ".shape",
            ".slider",
            ".smith",
            ".template",
            ".ternary",
            ".title",
            ".updatemenu",
            ".xaxis",
            ".yaxis",
        ],
        [
            "._activeselection.Activeselection",
            "._activeshape.Activeshape",
            "._annotation.Annotation",
            "._coloraxis.Coloraxis",
            "._colorscale.Colorscale",
            "._font.Font",
            "._geo.Geo",
            "._grid.Grid",
            "._hoverlabel.Hoverlabel",
            "._image.Image",
            "._legend.Legend",
            "._map.Map",
            "._mapbox.Mapbox",
            "._margin.Margin",
            "._modebar.Modebar",
            "._newselection.Newselection",
            "._newshape.Newshape",
            "._polar.Polar",
            "._scene.Scene",
            "._selection.Selection",
            "._shape.Shape",
            "._slider.Slider",
            "._smith.Smith",
            "._template.Template",
            "._ternary.Ternary",
            "._title.Title",
            "._transition.Transition",
            "._uniformtext.Uniformtext",
            "._updatemenu.Updatemenu",
            "._xaxis.XAxis",
            "._yaxis.YAxis",
        ],
    )
