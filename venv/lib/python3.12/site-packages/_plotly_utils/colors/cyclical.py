"""
Cyclical color scales are appropriate for continuous data that has a natural cyclical \
structure, such as temporal data (hour of day, day of week, day of year, seasons) or
complex numbers or other phase data.
"""

from ._swatches import _swatches, _swatches_continuous, _swatches_cyclical


def swatches(template=None):
    return _swatches(__name__, globals(), template)


swatches.__doc__ = _swatches.__doc__


def swatches_continuous(template=None):
    return _swatches_continuous(__name__, globals(), template)


swatches_continuous.__doc__ = _swatches_continuous.__doc__


def swatches_cyclical(template=None):
    return _swatches_cyclical(__name__, globals(), template)


swatches_cyclical.__doc__ = _swatches_cyclical.__doc__


Twilight = [
    "#e2d9e2",
    "#9ebbc9",
    "#6785be",
    "#5e43a5",
    "#421257",
    "#471340",
    "#8e2c50",
    "#ba6657",
    "#ceac94",
    "#e2d9e2",
]
IceFire = [
    "#000000",
    "#001f4d",
    "#003786",
    "#0e58a8",
    "#217eb8",
    "#30a4ca",
    "#54c8df",
    "#9be4ef",
    "#e1e9d1",
    "#f3d573",
    "#e7b000",
    "#da8200",
    "#c65400",
    "#ac2301",
    "#820000",
    "#4c0000",
    "#000000",
]
Edge = [
    "#313131",
    "#3d019d",
    "#3810dc",
    "#2d47f9",
    "#2593ff",
    "#2adef6",
    "#60fdfa",
    "#aefdff",
    "#f3f3f1",
    "#fffda9",
    "#fafd5b",
    "#f7da29",
    "#ff8e25",
    "#f8432d",
    "#d90d39",
    "#97023d",
    "#313131",
]
Phase = [
    "rgb(167, 119, 12)",
    "rgb(197, 96, 51)",
    "rgb(217, 67, 96)",
    "rgb(221, 38, 163)",
    "rgb(196, 59, 224)",
    "rgb(153, 97, 244)",
    "rgb(95, 127, 228)",
    "rgb(40, 144, 183)",
    "rgb(15, 151, 136)",
    "rgb(39, 153, 79)",
    "rgb(119, 141, 17)",
    "rgb(167, 119, 12)",
]
HSV = [
    "#ff0000",
    "#ffa700",
    "#afff00",
    "#08ff00",
    "#00ff9f",
    "#00b7ff",
    "#0010ff",
    "#9700ff",
    "#ff00bf",
    "#ff0000",
]
mrybm = [
    "#f884f7",
    "#f968c4",
    "#ea4388",
    "#cf244b",
    "#b51a15",
    "#bd4304",
    "#cc6904",
    "#d58f04",
    "#cfaa27",
    "#a19f62",
    "#588a93",
    "#2269c4",
    "#3e3ef0",
    "#6b4ef9",
    "#956bfa",
    "#cd7dfe",
    "#f884f7",
]
mygbm = [
    "#ef55f1",
    "#fb84ce",
    "#fbafa1",
    "#fcd471",
    "#f0ed35",
    "#c6e516",
    "#96d310",
    "#61c10b",
    "#31ac28",
    "#439064",
    "#3d719a",
    "#284ec8",
    "#2e21ea",
    "#6324f5",
    "#9139fa",
    "#c543fa",
    "#ef55f1",
]

Edge_r = Edge[::-1]
HSV_r = HSV[::-1]
IceFire_r = IceFire[::-1]
Phase_r = Phase[::-1]
Twilight_r = Twilight[::-1]
mrybm_r = mrybm[::-1]
mygbm_r = mygbm[::-1]

__all__ = [
    "swatches",
    "swatches_cyclical",
]
