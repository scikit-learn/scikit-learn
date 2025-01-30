import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._yaxis import YaxisValidator
    from ._xaxis import XaxisValidator
    from ._width import WidthValidator
    from ._waterfallmode import WaterfallmodeValidator
    from ._waterfallgroupgap import WaterfallgroupgapValidator
    from ._waterfallgap import WaterfallgapValidator
    from ._violinmode import ViolinmodeValidator
    from ._violingroupgap import ViolingroupgapValidator
    from ._violingap import ViolingapValidator
    from ._updatemenudefaults import UpdatemenudefaultsValidator
    from ._updatemenus import UpdatemenusValidator
    from ._uniformtext import UniformtextValidator
    from ._uirevision import UirevisionValidator
    from ._treemapcolorway import TreemapcolorwayValidator
    from ._transition import TransitionValidator
    from ._title import TitleValidator
    from ._ternary import TernaryValidator
    from ._template import TemplateValidator
    from ._sunburstcolorway import SunburstcolorwayValidator
    from ._spikedistance import SpikedistanceValidator
    from ._smith import SmithValidator
    from ._sliderdefaults import SliderdefaultsValidator
    from ._sliders import SlidersValidator
    from ._showlegend import ShowlegendValidator
    from ._shapedefaults import ShapedefaultsValidator
    from ._shapes import ShapesValidator
    from ._separators import SeparatorsValidator
    from ._selectiondefaults import SelectiondefaultsValidator
    from ._selections import SelectionsValidator
    from ._selectionrevision import SelectionrevisionValidator
    from ._selectdirection import SelectdirectionValidator
    from ._scene import SceneValidator
    from ._scattermode import ScattermodeValidator
    from ._scattergap import ScattergapValidator
    from ._polar import PolarValidator
    from ._plot_bgcolor import Plot_BgcolorValidator
    from ._piecolorway import PiecolorwayValidator
    from ._paper_bgcolor import Paper_BgcolorValidator
    from ._newshape import NewshapeValidator
    from ._newselection import NewselectionValidator
    from ._modebar import ModebarValidator
    from ._minreducedwidth import MinreducedwidthValidator
    from ._minreducedheight import MinreducedheightValidator
    from ._metasrc import MetasrcValidator
    from ._meta import MetaValidator
    from ._margin import MarginValidator
    from ._mapbox import MapboxValidator
    from ._map import MapValidator
    from ._legend import LegendValidator
    from ._imagedefaults import ImagedefaultsValidator
    from ._images import ImagesValidator
    from ._iciclecolorway import IciclecolorwayValidator
    from ._hoversubplots import HoversubplotsValidator
    from ._hovermode import HovermodeValidator
    from ._hoverlabel import HoverlabelValidator
    from ._hoverdistance import HoverdistanceValidator
    from ._hidesources import HidesourcesValidator
    from ._hiddenlabelssrc import HiddenlabelssrcValidator
    from ._hiddenlabels import HiddenlabelsValidator
    from ._height import HeightValidator
    from ._grid import GridValidator
    from ._geo import GeoValidator
    from ._funnelmode import FunnelmodeValidator
    from ._funnelgroupgap import FunnelgroupgapValidator
    from ._funnelgap import FunnelgapValidator
    from ._funnelareacolorway import FunnelareacolorwayValidator
    from ._font import FontValidator
    from ._extendtreemapcolors import ExtendtreemapcolorsValidator
    from ._extendsunburstcolors import ExtendsunburstcolorsValidator
    from ._extendpiecolors import ExtendpiecolorsValidator
    from ._extendiciclecolors import ExtendiciclecolorsValidator
    from ._extendfunnelareacolors import ExtendfunnelareacolorsValidator
    from ._editrevision import EditrevisionValidator
    from ._dragmode import DragmodeValidator
    from ._datarevision import DatarevisionValidator
    from ._computed import ComputedValidator
    from ._colorway import ColorwayValidator
    from ._colorscale import ColorscaleValidator
    from ._coloraxis import ColoraxisValidator
    from ._clickmode import ClickmodeValidator
    from ._calendar import CalendarValidator
    from ._boxmode import BoxmodeValidator
    from ._boxgroupgap import BoxgroupgapValidator
    from ._boxgap import BoxgapValidator
    from ._barnorm import BarnormValidator
    from ._barmode import BarmodeValidator
    from ._bargroupgap import BargroupgapValidator
    from ._bargap import BargapValidator
    from ._barcornerradius import BarcornerradiusValidator
    from ._autotypenumbers import AutotypenumbersValidator
    from ._autosize import AutosizeValidator
    from ._annotationdefaults import AnnotationdefaultsValidator
    from ._annotations import AnnotationsValidator
    from ._activeshape import ActiveshapeValidator
    from ._activeselection import ActiveselectionValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._yaxis.YaxisValidator",
            "._xaxis.XaxisValidator",
            "._width.WidthValidator",
            "._waterfallmode.WaterfallmodeValidator",
            "._waterfallgroupgap.WaterfallgroupgapValidator",
            "._waterfallgap.WaterfallgapValidator",
            "._violinmode.ViolinmodeValidator",
            "._violingroupgap.ViolingroupgapValidator",
            "._violingap.ViolingapValidator",
            "._updatemenudefaults.UpdatemenudefaultsValidator",
            "._updatemenus.UpdatemenusValidator",
            "._uniformtext.UniformtextValidator",
            "._uirevision.UirevisionValidator",
            "._treemapcolorway.TreemapcolorwayValidator",
            "._transition.TransitionValidator",
            "._title.TitleValidator",
            "._ternary.TernaryValidator",
            "._template.TemplateValidator",
            "._sunburstcolorway.SunburstcolorwayValidator",
            "._spikedistance.SpikedistanceValidator",
            "._smith.SmithValidator",
            "._sliderdefaults.SliderdefaultsValidator",
            "._sliders.SlidersValidator",
            "._showlegend.ShowlegendValidator",
            "._shapedefaults.ShapedefaultsValidator",
            "._shapes.ShapesValidator",
            "._separators.SeparatorsValidator",
            "._selectiondefaults.SelectiondefaultsValidator",
            "._selections.SelectionsValidator",
            "._selectionrevision.SelectionrevisionValidator",
            "._selectdirection.SelectdirectionValidator",
            "._scene.SceneValidator",
            "._scattermode.ScattermodeValidator",
            "._scattergap.ScattergapValidator",
            "._polar.PolarValidator",
            "._plot_bgcolor.Plot_BgcolorValidator",
            "._piecolorway.PiecolorwayValidator",
            "._paper_bgcolor.Paper_BgcolorValidator",
            "._newshape.NewshapeValidator",
            "._newselection.NewselectionValidator",
            "._modebar.ModebarValidator",
            "._minreducedwidth.MinreducedwidthValidator",
            "._minreducedheight.MinreducedheightValidator",
            "._metasrc.MetasrcValidator",
            "._meta.MetaValidator",
            "._margin.MarginValidator",
            "._mapbox.MapboxValidator",
            "._map.MapValidator",
            "._legend.LegendValidator",
            "._imagedefaults.ImagedefaultsValidator",
            "._images.ImagesValidator",
            "._iciclecolorway.IciclecolorwayValidator",
            "._hoversubplots.HoversubplotsValidator",
            "._hovermode.HovermodeValidator",
            "._hoverlabel.HoverlabelValidator",
            "._hoverdistance.HoverdistanceValidator",
            "._hidesources.HidesourcesValidator",
            "._hiddenlabelssrc.HiddenlabelssrcValidator",
            "._hiddenlabels.HiddenlabelsValidator",
            "._height.HeightValidator",
            "._grid.GridValidator",
            "._geo.GeoValidator",
            "._funnelmode.FunnelmodeValidator",
            "._funnelgroupgap.FunnelgroupgapValidator",
            "._funnelgap.FunnelgapValidator",
            "._funnelareacolorway.FunnelareacolorwayValidator",
            "._font.FontValidator",
            "._extendtreemapcolors.ExtendtreemapcolorsValidator",
            "._extendsunburstcolors.ExtendsunburstcolorsValidator",
            "._extendpiecolors.ExtendpiecolorsValidator",
            "._extendiciclecolors.ExtendiciclecolorsValidator",
            "._extendfunnelareacolors.ExtendfunnelareacolorsValidator",
            "._editrevision.EditrevisionValidator",
            "._dragmode.DragmodeValidator",
            "._datarevision.DatarevisionValidator",
            "._computed.ComputedValidator",
            "._colorway.ColorwayValidator",
            "._colorscale.ColorscaleValidator",
            "._coloraxis.ColoraxisValidator",
            "._clickmode.ClickmodeValidator",
            "._calendar.CalendarValidator",
            "._boxmode.BoxmodeValidator",
            "._boxgroupgap.BoxgroupgapValidator",
            "._boxgap.BoxgapValidator",
            "._barnorm.BarnormValidator",
            "._barmode.BarmodeValidator",
            "._bargroupgap.BargroupgapValidator",
            "._bargap.BargapValidator",
            "._barcornerradius.BarcornerradiusValidator",
            "._autotypenumbers.AutotypenumbersValidator",
            "._autosize.AutosizeValidator",
            "._annotationdefaults.AnnotationdefaultsValidator",
            "._annotations.AnnotationsValidator",
            "._activeshape.ActiveshapeValidator",
            "._activeselection.ActiveselectionValidator",
        ],
    )
