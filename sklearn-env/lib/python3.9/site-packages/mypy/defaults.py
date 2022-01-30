import os

from typing_extensions import Final

PYTHON2_VERSION: Final = (2, 7)
PYTHON3_VERSION: Final = (3, 6)
PYTHON3_VERSION_MIN: Final = (3, 4)
CACHE_DIR: Final = ".mypy_cache"
CONFIG_FILE: Final = ["mypy.ini", ".mypy.ini"]
PYPROJECT_CONFIG_FILES: Final = [
    "pyproject.toml",
]
SHARED_CONFIG_FILES: Final = [
    "setup.cfg",
]
USER_CONFIG_FILES: Final = [
    "~/.config/mypy/config",
    "~/.mypy.ini",
]
if os.environ.get("XDG_CONFIG_HOME"):
    USER_CONFIG_FILES.insert(0, os.path.join(os.environ["XDG_CONFIG_HOME"], "mypy/config"))

CONFIG_FILES: Final = (
    CONFIG_FILE + PYPROJECT_CONFIG_FILES + SHARED_CONFIG_FILES + USER_CONFIG_FILES
)

# This must include all reporters defined in mypy.report. This is defined here
# to make reporter names available without importing mypy.report -- this speeds
# up startup.
REPORTER_NAMES: Final = [
    "linecount",
    "any-exprs",
    "linecoverage",
    "memory-xml",
    "cobertura-xml",
    "xml",
    "xslt-html",
    "xslt-txt",
    "html",
    "txt",
    "lineprecision",
]

# Threshold after which we sometimes filter out most errors to avoid very
# verbose output
MANY_ERRORS_THRESHOLD: Final = 200
