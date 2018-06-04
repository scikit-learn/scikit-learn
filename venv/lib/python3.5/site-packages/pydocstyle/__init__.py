from .checker import check
from .violations import Error, conventions
from .utils import __version__

# Temporary hotfix for flake8-docstrings
from .checker import ConventionChecker, tokenize_open
from .parser import AllError
