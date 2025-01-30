"""
Diverging color scales are appropriate for continuous data that has a natural midpoint \
other otherwise informative special value, such as 0 altitude, or the boiling point
of a liquid. The color scales in this module are \
mostly meant to be passed in as the `color_continuous_scale` argument to various \
functions, and to be used with the `color_continuous_midpoint` argument.
"""

from .colorbrewer import (  # noqa: F401
    BrBG,
    PRGn,
    PiYG,
    PuOr,
    RdBu,
    RdGy,
    RdYlBu,
    RdYlGn,
    Spectral,
    BrBG_r,
    PRGn_r,
    PiYG_r,
    PuOr_r,
    RdBu_r,
    RdGy_r,
    RdYlBu_r,
    RdYlGn_r,
    Spectral_r,
)
from .cmocean import (  # noqa: F401
    balance,
    delta,
    curl,
    oxy,
    balance_r,
    delta_r,
    curl_r,
    oxy_r,
)
from .carto import (  # noqa: F401
    Armyrose,
    Fall,
    Geyser,
    Temps,
    Tealrose,
    Tropic,
    Earth,
    Armyrose_r,
    Fall_r,
    Geyser_r,
    Temps_r,
    Tealrose_r,
    Tropic_r,
    Earth_r,
)

from .plotlyjs import Picnic, Portland, Picnic_r, Portland_r  # noqa: F401

from ._swatches import _swatches, _swatches_continuous


def swatches(template=None):
    return _swatches(__name__, globals(), template)


swatches.__doc__ = _swatches.__doc__


def swatches_continuous(template=None):
    return _swatches_continuous(__name__, globals(), template)


swatches_continuous.__doc__ = _swatches_continuous.__doc__


__all__ = ["swatches"]
