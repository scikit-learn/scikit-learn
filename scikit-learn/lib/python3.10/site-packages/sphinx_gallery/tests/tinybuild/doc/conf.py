"""Sphinx configuration for tinybuild."""

import os.path as op

import sphinx


class CustomSortKey:
    """Basically a clone of FilenameSortKey, for testing custom class sorters."""

    def __init__(self, src_dir):
        self.src_dir = src_dir

    def __call__(self, item):
        """Provide the custom sort key."""
        return item


# Where our helpers live
util_root = "sphinx_gallery.tests.tinybuild.utils"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx_gallery.gen_gallery",
    "sphinx.ext.graphviz",
    "jupyterlite_sphinx",
]
templates_path = ["_templates"]
autosummary_generate = True
source_suffix = ".rst"
master_doc = "index"
exclude_patterns = ["_build"]
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
    "joblib": ("https://joblib.readthedocs.io/en/latest", None),
}
html_sidebars = {
    "**": [
        "globaltoc.html",
        "searchbox.html",
        "sg_download_links.html",
        "sg_launcher_links.html",
    ]
}

sphinx_gallery_conf = {
    "doc_module": ("sphinx_gallery",),
    "prefer_full_module": {r"sphinx_gallery\._dummy"},
    "reference_url": {
        "sphinx_gallery": None,
        "scipy": "http://docs.scipy.org/doc/scipy/wrong_url",  # bad one
    },
    "binder": {
        "org": "sphinx-gallery",
        "repo": "sphinx-gallery.github.io",
        "branch": "master",
        "binderhub_url": "https://mybinder.org",
        "dependencies": "./binder/requirements.txt",
        "notebooks_dir": "notebooks",
        "use_jupyter_lab": True,
    },
    "jupyterlite": {
        "notebook_modification_function": f"{util_root}.notebook_modification_function",
    },
    "examples_dirs": [
        "../examples/",
        "../examples_with_rst/",
        "../examples_rst_index",
        "../examples_README_header",
    ],
    "example_extensions": {".py", ".cpp", ".m", ".jl"},
    "filetype_parsers": {".m": "Matlab"},
    "reset_argv": f"{util_root}.reset_argv",
    "reset_modules": (f"{util_root}.mock_scrape_problem", "matplotlib"),
    "gallery_dirs": [
        "auto_examples",
        "auto_examples_with_rst",
        "auto_examples_rst_index",
        "auto_examples_README_header",
    ],
    "backreferences_dir": "gen_modules/backreferences",
    "subsection_order": f"{util_root}.noop_key",
    "within_subsection_order": CustomSortKey,
    "minigallery_sort_order": f"{util_root}.mock_sort",
    "image_scrapers": (f"{util_root}.matplotlib_format_scraper",),
    "expected_failing_examples": [
        "../examples/future/plot_future_imports_broken.py",
        "../examples/plot_scraper_broken.py",
        "../examples/plot_failing_example.py",
        "../examples/plot_failing_example_thumbnail.py",
    ],
    "show_memory": False,
    "compress_images": ("images", "thumbnails"),
    "junit": op.join("sphinx-gallery", "junit-results.xml"),
    "matplotlib_animations": (True, "mp4"),
    "pypandoc": True,
    "image_srcset": ["2x"],
    "exclude_implicit_doc": ["figure_rst"],
    "show_api_usage": True,
    "copyfile_regex": r".*\.rst",
    "recommender": {"enable": True, "n_examples": 3},
    "parallel": 2,
}
nitpicky = True
highlight_language = "python3"
html_static_path = ["_static_nonstandard"]


class ReportChanged:
    """For logging files changed by connecting to Sphinx `env-get-outdated` event."""

    def __repr__(self):
        return "ReportChanged"

    def __call__(self, app, env, added, changed, removed):
        """Log files that have changed."""
        from sphinx.util.console import bold

        logger = sphinx.util.logging.getLogger("sphinx-gallery")
        if changed:
            logger.info(bold(f"\nFiles changed ({len(changed)}):"))
            for fname in sorted(changed):
                logger.info(f"     - {fname}")
        return []


def setup(app):
    """Setup Sphinx."""
    app.connect("env-get-outdated", ReportChanged())
