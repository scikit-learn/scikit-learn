# This file is not meant for public use and will be removed in SciPy v2.0.0.
# Use the `scipy.interpolate` namespace for importing the functions
# included below.

from scipy._lib.deprecation import _sub_module_deprecation


__all__ = [  # noqa: F822
    'bispeu',
    'bispev',
    'curfit',
    'dblint',
    'fpchec',
    'fpcurf0',
    'fpcurf1',
    'fpcurfm1',
    'parcur',
    'parder',
    'pardeu',
    'pardtc',
    'percur',
    'regrid_smth',
    'regrid_smth_spher',
    'spalde',
    'spherfit_lsq',
    'spherfit_smth',
    'splder',
    'splev',
    'splint',
    'sproot',
    'surfit_lsq',
    'surfit_smth',
    'types',
]


def __dir__():
    return __all__


def __getattr__(name):
    return _sub_module_deprecation(sub_package="interpolate", module="dfitpack",
                                   private_modules=["_dfitpack"], all=__all__,
                                   attribute=name)
