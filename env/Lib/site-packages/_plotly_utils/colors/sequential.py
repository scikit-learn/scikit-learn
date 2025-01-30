"""
Sequential color scales are appropriate for most continuous data, but in some cases it \
can be helpful to use a `plotly.colors.diverging` or \
`plotly.colors.cyclical` scale instead. The color scales in this module are \
mostly meant to be passed in as the `color_continuous_scale` argument to various functions.
"""

from ._swatches import _swatches, _swatches_continuous


def swatches(template=None):
    return _swatches(__name__, globals(), template)


swatches.__doc__ = _swatches.__doc__


def swatches_continuous(template=None):
    return _swatches_continuous(__name__, globals(), template)


swatches_continuous.__doc__ = _swatches_continuous.__doc__

Plotly3 = [
    "#0508b8",
    "#1910d8",
    "#3c19f0",
    "#6b1cfb",
    "#981cfd",
    "#bf1cfd",
    "#dd2bfd",
    "#f246fe",
    "#fc67fd",
    "#fe88fc",
    "#fea5fd",
    "#febefe",
    "#fec3fe",
]

Viridis = [
    "#440154",
    "#482878",
    "#3e4989",
    "#31688e",
    "#26828e",
    "#1f9e89",
    "#35b779",
    "#6ece58",
    "#b5de2b",
    "#fde725",
]
Cividis = [
    "#00224e",
    "#123570",
    "#3b496c",
    "#575d6d",
    "#707173",
    "#8a8678",
    "#a59c74",
    "#c3b369",
    "#e1cc55",
    "#fee838",
]

Inferno = [
    "#000004",
    "#1b0c41",
    "#4a0c6b",
    "#781c6d",
    "#a52c60",
    "#cf4446",
    "#ed6925",
    "#fb9b06",
    "#f7d13d",
    "#fcffa4",
]
Magma = [
    "#000004",
    "#180f3d",
    "#440f76",
    "#721f81",
    "#9e2f7f",
    "#cd4071",
    "#f1605d",
    "#fd9668",
    "#feca8d",
    "#fcfdbf",
]
Plasma = [
    "#0d0887",
    "#46039f",
    "#7201a8",
    "#9c179e",
    "#bd3786",
    "#d8576b",
    "#ed7953",
    "#fb9f3a",
    "#fdca26",
    "#f0f921",
]
Turbo = [
    "#30123b",
    "#4145ab",
    "#4675ed",
    "#39a2fc",
    "#1bcfd4",
    "#24eca6",
    "#61fc6c",
    "#a4fc3b",
    "#d1e834",
    "#f3c63a",
    "#fe9b2d",
    "#f36315",
    "#d93806",
    "#b11901",
    "#7a0402",
]

Cividis_r = Cividis[::-1]
Inferno_r = Inferno[::-1]
Magma_r = Magma[::-1]
Plasma_r = Plasma[::-1]
Plotly3_r = Plotly3[::-1]
Turbo_r = Turbo[::-1]
Viridis_r = Viridis[::-1]

from .plotlyjs import (  # noqa: F401
    Blackbody,
    Bluered,
    Electric,
    Hot,
    Jet,
    Rainbow,
    Blackbody_r,
    Bluered_r,
    Electric_r,
    Hot_r,
    Jet_r,
    Rainbow_r,
)

from .colorbrewer import (  # noqa: F401
    Blues,
    BuGn,
    BuPu,
    GnBu,
    Greens,
    Greys,
    OrRd,
    Oranges,
    PuBu,
    PuBuGn,
    PuRd,
    Purples,
    RdBu,
    RdPu,
    Reds,
    YlGn,
    YlGnBu,
    YlOrBr,
    YlOrRd,
    Blues_r,
    BuGn_r,
    BuPu_r,
    GnBu_r,
    Greens_r,
    Greys_r,
    OrRd_r,
    Oranges_r,
    PuBu_r,
    PuBuGn_r,
    PuRd_r,
    Purples_r,
    RdBu_r,
    RdPu_r,
    Reds_r,
    YlGn_r,
    YlGnBu_r,
    YlOrBr_r,
    YlOrRd_r,
)

from .cmocean import (  # noqa: F401
    turbid,
    thermal,
    haline,
    solar,
    ice,
    gray,
    deep,
    dense,
    algae,
    matter,
    speed,
    amp,
    tempo,
    turbid_r,
    thermal_r,
    haline_r,
    solar_r,
    ice_r,
    gray_r,
    deep_r,
    dense_r,
    algae_r,
    matter_r,
    speed_r,
    amp_r,
    tempo_r,
)

from .carto import (  # noqa: F401
    Burg,
    Burgyl,
    Redor,
    Oryel,
    Peach,
    Pinkyl,
    Mint,
    Blugrn,
    Darkmint,
    Emrld,
    Aggrnyl,
    Bluyl,
    Teal,
    Tealgrn,
    Purp,
    Purpor,
    Sunset,
    Magenta,
    Sunsetdark,
    Agsunset,
    Brwnyl,
    Burg_r,
    Burgyl_r,
    Redor_r,
    Oryel_r,
    Peach_r,
    Pinkyl_r,
    Mint_r,
    Blugrn_r,
    Darkmint_r,
    Emrld_r,
    Aggrnyl_r,
    Bluyl_r,
    Teal_r,
    Tealgrn_r,
    Purp_r,
    Purpor_r,
    Sunset_r,
    Magenta_r,
    Sunsetdark_r,
    Agsunset_r,
    Brwnyl_r,
)

__all__ = ["swatches"]
