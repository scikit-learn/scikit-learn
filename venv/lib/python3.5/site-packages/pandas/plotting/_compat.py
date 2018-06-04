# being a bit too dynamic
# pylint: disable=E1101
from __future__ import division
import operator

from distutils.version import LooseVersion


def _mpl_version(version, op):
    def inner():
        try:
            import matplotlib as mpl
        except ImportError:
            return False
        return (op(LooseVersion(mpl.__version__), LooseVersion(version)) and
                str(mpl.__version__)[0] != '0')

    return inner


_mpl_ge_1_2_1 = _mpl_version('1.2.1', operator.ge)
_mpl_le_1_2_1 = _mpl_version('1.2.1', operator.le)
_mpl_ge_1_3_1 = _mpl_version('1.3.1', operator.ge)
_mpl_ge_1_4_0 = _mpl_version('1.4.0', operator.ge)
_mpl_ge_1_4_1 = _mpl_version('1.4.1', operator.ge)
_mpl_ge_1_5_0 = _mpl_version('1.5.0', operator.ge)
_mpl_ge_2_0_0 = _mpl_version('2.0.0', operator.ge)
_mpl_le_2_0_0 = _mpl_version('2.0.0', operator.le)
_mpl_ge_2_0_1 = _mpl_version('2.0.1', operator.ge)
_mpl_ge_2_1_0 = _mpl_version('2.1.0', operator.ge)
_mpl_ge_2_2_0 = _mpl_version('2.2.0', operator.ge)
