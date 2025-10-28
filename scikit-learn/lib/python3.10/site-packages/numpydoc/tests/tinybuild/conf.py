import os
import sys

path = os.path.dirname(__file__)
if path not in sys.path:
    sys.path.insert(0, path)
import numpydoc_test_module

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "numpydoc",
]
project = "numpydoc_test_module"
autosummary_generate = True
autodoc_default_options = {"inherited-members": None}
source_suffix = ".rst"
exclude_patterns = ["_build"]
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
}
nitpicky = True
highlight_language = "python3"
numpydoc_class_members_toctree = False
numpydoc_xref_param_type = True
