"""Callables referenced by ``sphinx_gallery_conf``.

Sphinx-Gallery resolves these by fully qualified name so that no function or class
object ends up in ``conf.py``'s config values, which would make the config
unpicklable and thereby block parallel reads. See
https://sphinx-gallery.github.io/stable/configuration.html#importing-callables
"""

import os
import re

from sphinx_gallery.notebook import add_code_cell, add_markdown_cell
from sphinx_gallery.sorting import ExampleTitleSortKey

import sklearn
from sklearn.externals._packaging.version import parse

EXAMPLES_DIR = "../examples"

_subsection_title = re.compile(r"^([\w ]+)\n-", re.MULTILINE)

_parsed_version = parse(sklearn.__version__)
if _parsed_version.is_postrelease:
    _release = _parsed_version.base_version
else:
    _release = sklearn.__version__

_default_global_config = sklearn.get_config()


def subsection_order(directory):
    """Sort example gallery by title of subsection.

    Assumes README.txt exists for all subsections and uses the subsection with
    dashes, '---', as the adornment.
    """
    src_path = os.path.normpath(os.path.join(EXAMPLES_DIR, directory))

    # Forces Release Highlights to the top
    if os.path.basename(src_path) == "release_highlights":
        return "0"

    readme = os.path.join(src_path, "README.txt")

    try:
        with open(readme, "r") as f:
            content = f.read()
    except FileNotFoundError:
        return directory

    title_match = _subsection_title.search(content)
    if title_match is not None:
        return title_match.group(1)
    return directory


class SKExampleTitleSortKey(ExampleTitleSortKey):
    """Sorts release highlights based on version number."""

    def __call__(self, filename):
        title = super().__call__(filename)
        prefix = "plot_release_highlights_"

        # Use title to sort if not a release highlight
        if not str(filename).startswith(prefix):
            return title

        major_minor = filename[len(prefix) :].split("_")[:2]
        version_float = float(".".join(major_minor))

        # negate to place the newest version highlights first
        return -version_float


def notebook_modification_function(notebook_content, notebook_filename):
    notebook_content_str = str(notebook_content)
    warning_template = "\n".join(
        [
            "<div class='alert alert-{message_class}'>",
            "",
            "# JupyterLite warning",
            "",
            "{message}",
            "</div>",
        ]
    )

    message_class = "warning"
    message = (
        "Running the scikit-learn examples in JupyterLite is experimental and you may"
        " encounter some unexpected behavior.\n\nThe main difference is that imports"
        " will take a lot longer than usual, for example the first `import sklearn` can"
        " take roughly 10-20s.\n\nIf you notice problems, feel free to open an"
        " [issue](https://github.com/scikit-learn/scikit-learn/issues/new/choose)"
        " about it."
    )

    markdown = warning_template.format(message_class=message_class, message=message)

    dummy_notebook_content = {"cells": []}
    add_markdown_cell(dummy_notebook_content, markdown)

    code_lines = []

    if "seaborn" in notebook_content_str:
        code_lines.append("%pip install seaborn")
    if "plotly.express" in notebook_content_str:
        code_lines.append("%pip install plotly nbformat")
    if "skimage" in notebook_content_str:
        code_lines.append("%pip install scikit-image")
    if "polars" in notebook_content_str:
        code_lines.append("%pip install polars")
    if "fetch_" in notebook_content_str:
        code_lines.extend(
            [
                "%pip install pyodide-http",
                "import pyodide_http",
                "pyodide_http.patch_all()",
            ]
        )
    # always import matplotlib and pandas to avoid Pyodide limitation with
    # imports inside functions
    code_lines.extend(["import matplotlib", "import pandas"])

    # Work around https://github.com/jupyterlite/pyodide-kernel/issues/166
    # and https://github.com/pyodide/micropip/issues/223 by installing the
    # dependencies first, and then scikit-learn from Anaconda.org.
    if "dev" in _release:
        dev_docs_specific_code = [
            "import piplite",
            "import joblib",
            "import threadpoolctl",
            "import scipy",
            "await piplite.install(\n"
            f"  'scikit-learn=={_release}',\n"
            "   index_urls='https://pypi.anaconda.org/scientific-python-nightly-wheels/simple',\n"
            ")",
        ]

        code_lines.extend(dev_docs_specific_code)

    if code_lines:
        code_lines = ["# JupyterLite-specific code"] + code_lines
        code = "\n".join(code_lines)
        add_code_cell(dummy_notebook_content, code)

    notebook_content["cells"] = (
        dummy_notebook_content["cells"] + notebook_content["cells"]
    )


def reset_sklearn_config(gallery_conf, fname):
    """Reset sklearn config to default values."""
    sklearn.set_config(**_default_global_config)
