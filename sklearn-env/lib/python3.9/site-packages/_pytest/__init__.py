__all__ = ["__version__"]

try:
    from ._version import version as __version__
except ImportError:
    # broken installation, we don't even try
    # unknown only works because we do poor mans version compare
    __version__ = "unknown"
