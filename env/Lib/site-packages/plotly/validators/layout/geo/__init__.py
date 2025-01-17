import sys
from typing import TYPE_CHECKING

if sys.version_info < (3, 7) or TYPE_CHECKING:
    from ._visible import VisibleValidator
    from ._uirevision import UirevisionValidator
    from ._subunitwidth import SubunitwidthValidator
    from ._subunitcolor import SubunitcolorValidator
    from ._showsubunits import ShowsubunitsValidator
    from ._showrivers import ShowriversValidator
    from ._showocean import ShowoceanValidator
    from ._showland import ShowlandValidator
    from ._showlakes import ShowlakesValidator
    from ._showframe import ShowframeValidator
    from ._showcountries import ShowcountriesValidator
    from ._showcoastlines import ShowcoastlinesValidator
    from ._scope import ScopeValidator
    from ._riverwidth import RiverwidthValidator
    from ._rivercolor import RivercolorValidator
    from ._resolution import ResolutionValidator
    from ._projection import ProjectionValidator
    from ._oceancolor import OceancolorValidator
    from ._lonaxis import LonaxisValidator
    from ._lataxis import LataxisValidator
    from ._landcolor import LandcolorValidator
    from ._lakecolor import LakecolorValidator
    from ._framewidth import FramewidthValidator
    from ._framecolor import FramecolorValidator
    from ._fitbounds import FitboundsValidator
    from ._domain import DomainValidator
    from ._countrywidth import CountrywidthValidator
    from ._countrycolor import CountrycolorValidator
    from ._coastlinewidth import CoastlinewidthValidator
    from ._coastlinecolor import CoastlinecolorValidator
    from ._center import CenterValidator
    from ._bgcolor import BgcolorValidator
else:
    from _plotly_utils.importers import relative_import

    __all__, __getattr__, __dir__ = relative_import(
        __name__,
        [],
        [
            "._visible.VisibleValidator",
            "._uirevision.UirevisionValidator",
            "._subunitwidth.SubunitwidthValidator",
            "._subunitcolor.SubunitcolorValidator",
            "._showsubunits.ShowsubunitsValidator",
            "._showrivers.ShowriversValidator",
            "._showocean.ShowoceanValidator",
            "._showland.ShowlandValidator",
            "._showlakes.ShowlakesValidator",
            "._showframe.ShowframeValidator",
            "._showcountries.ShowcountriesValidator",
            "._showcoastlines.ShowcoastlinesValidator",
            "._scope.ScopeValidator",
            "._riverwidth.RiverwidthValidator",
            "._rivercolor.RivercolorValidator",
            "._resolution.ResolutionValidator",
            "._projection.ProjectionValidator",
            "._oceancolor.OceancolorValidator",
            "._lonaxis.LonaxisValidator",
            "._lataxis.LataxisValidator",
            "._landcolor.LandcolorValidator",
            "._lakecolor.LakecolorValidator",
            "._framewidth.FramewidthValidator",
            "._framecolor.FramecolorValidator",
            "._fitbounds.FitboundsValidator",
            "._domain.DomainValidator",
            "._countrywidth.CountrywidthValidator",
            "._countrycolor.CountrycolorValidator",
            "._coastlinewidth.CoastlinewidthValidator",
            "._coastlinecolor.CoastlinecolorValidator",
            "._center.CenterValidator",
            "._bgcolor.BgcolorValidator",
        ],
    )
