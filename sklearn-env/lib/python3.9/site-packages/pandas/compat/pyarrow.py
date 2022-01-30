""" support pyarrow compatibility across versions """

from pandas.util.version import Version

try:
    import pyarrow as pa

    _pa_version = pa.__version__
    _palv = Version(_pa_version)
    pa_version_under1p01 = _palv < Version("1.0.1")
    pa_version_under2p0 = _palv < Version("2.0.0")
    pa_version_under3p0 = _palv < Version("3.0.0")
    pa_version_under4p0 = _palv < Version("4.0.0")
    pa_version_under5p0 = _palv < Version("5.0.0")
    pa_version_under6p0 = _palv < Version("6.0.0")
except ImportError:
    pa_version_under1p01 = True
    pa_version_under2p0 = True
    pa_version_under3p0 = True
    pa_version_under4p0 = True
    pa_version_under5p0 = True
    pa_version_under6p0 = True
