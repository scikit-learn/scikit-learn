# scikit-learn documentation build configuration file, created by
# sphinx-quickstart on Fri Jan  8 09:13:42 2010.
#
# This file is execfile()d with the current directory set to its containing
# dir.
#
# Note that not all possible configuration values are present in this
# autogenerated file.
#
# All configuration values have a default; values that are commented out
# serve to show the default.

import json
import os
import re
import sys
import warnings
from datetime import datetime
from pathlib import Path
from urllib.request import urlopen

from sklearn.externals._packaging.version import parse
from sklearn.utils._testing import turn_warnings_into_errors

# If extensions (or modules to document with autodoc) are in another
# directory, add these directories to sys.path here. If the directory
# is relative to the documentation root, use os.path.abspath to make it
# absolute, like shown here.
sys.path.insert(0, os.path.abspath("."))
sys.path.insert(0, os.path.abspath("sphinxext"))

import jinja2
import sphinx_gallery
from github_link import make_linkcode_resolve
from sphinx.util.logging import getLogger
from sphinx_gallery.notebook import add_code_cell, add_markdown_cell
from sphinx_gallery.sorting import ExampleTitleSortKey

logger = getLogger(__name__)

try:
    # Configure plotly to integrate its output into the HTML pages generated by
    # sphinx-gallery.
    import plotly.io as pio

    pio.renderers.default = "sphinx_gallery"
except ImportError:
    # Make it possible to render the doc when not running the examples
    # that need plotly.
    pass

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "numpydoc",
    "sphinx.ext.linkcode",
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "sphinx.ext.imgconverter",
    "sphinx_gallery.gen_gallery",
    "sphinx-prompt",
    "sphinx_copybutton",
    "sphinxext.opengraph",
    "matplotlib.sphinxext.plot_directive",
    "sphinxcontrib.sass",
    "sphinx_remove_toctrees",
    "sphinx_design",
    # See sphinxext/
    "allow_nan_estimators",
    "autoshortsummary",
    "doi_role",
    "dropdown_anchors",
    "override_pst_pagetoc",
    "sphinx_issues",
]

# Specify how to identify the prompt when copying code snippets
copybutton_prompt_text = r">>> |\.\.\. "
copybutton_prompt_is_regexp = True
copybutton_exclude = "style"

try:
    import jupyterlite_sphinx  # noqa: F401

    extensions.append("jupyterlite_sphinx")
    with_jupyterlite = True
except ImportError:
    # In some cases we don't want to require jupyterlite_sphinx to be installed,
    # e.g. the doc-min-dependencies build
    warnings.warn(
        "jupyterlite_sphinx is not installed, you need to install it "
        "if you want JupyterLite links to appear in each example"
    )
    with_jupyterlite = False

# Produce `plot::` directives for examples that contain `import matplotlib` or
# `from matplotlib import`.
numpydoc_use_plots = True

# Options for the `::plot` directive:
# https://matplotlib.org/stable/api/sphinxext_plot_directive_api.html
plot_formats = ["png"]
plot_include_source = True
plot_html_show_formats = False
plot_html_show_source_link = False

# We do not need the table of class members because `sphinxext/override_pst_pagetoc.py`
# will show them in the secondary sidebar
numpydoc_show_class_members = False
numpydoc_show_inherited_class_members = False

# We want in-page toc of class members instead of a separate page for each entry
numpydoc_class_members_toctree = False


# For maths, use mathjax by default and svg if NO_MATHJAX env variable is set
# (useful for viewing the doc offline)
if os.environ.get("NO_MATHJAX"):
    extensions.append("sphinx.ext.imgmath")
    imgmath_image_format = "svg"
    mathjax_path = ""
else:
    extensions.append("sphinx.ext.mathjax")
    mathjax_path = "https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js"

# Add any paths that contain templates here, relative to this directory.
templates_path = ["templates"]

# generate autosummary even if no references
autosummary_generate = True

# The suffix of source filenames.
source_suffix = ".rst"

# The encoding of source files.
source_encoding = "utf-8"

# The main toctree document.
root_doc = "index"

# General information about the project.
project = "scikit-learn"
copyright = f"2007 - {datetime.now().year}, scikit-learn developers (BSD License)"

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
import sklearn

parsed_version = parse(sklearn.__version__)
version = ".".join(parsed_version.base_version.split(".")[:2])
# The full version, including alpha/beta/rc tags.
# Removes post from release name
if parsed_version.is_postrelease:
    release = parsed_version.base_version
else:
    release = sklearn.__version__

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
# language = None

# There are two options for replacing |today|: either, you set today to some
# non-false value, then it is used:
# today = ''
# Else, today_fmt is used as the format for a strftime call.
# today_fmt = '%B %d, %Y'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = [
    "_build",
    "templates",
    "includes",
    "**/sg_execution_times.rst",
]

# The reST default role (used for this markup: `text`) to use for all
# documents.
default_role = "literal"

# If true, '()' will be appended to :func: etc. cross-reference text.
add_function_parentheses = False

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
# add_module_names = True

# If true, sectionauthor and moduleauthor directives will be shown in the
# output. They are ignored by default.
# show_authors = False

# A list of ignored prefixes for module index sorting.
# modindex_common_prefix = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  Major themes that come with
# Sphinx are currently 'default' and 'sphinxdoc'.
html_theme = "pydata_sphinx_theme"

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
html_theme_options = {
    # -- General configuration ------------------------------------------------
    "sidebar_includehidden": True,
    "use_edit_page_button": True,
    "external_links": [],
    "icon_links_label": "Icon Links",
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/scikit-learn/scikit-learn",
            "icon": "fa-brands fa-square-github",
            "type": "fontawesome",
        },
    ],
    "analytics": {
        "plausible_analytics_domain": "scikit-learn.org",
        "plausible_analytics_url": "https://views.scientific-python.org/js/script.js",
    },
    # If "prev-next" is included in article_footer_items, then setting show_prev_next
    # to True would repeat prev and next links. See
    # https://github.com/pydata/pydata-sphinx-theme/blob/b731dc230bc26a3d1d1bb039c56c977a9b3d25d8/src/pydata_sphinx_theme/theme/pydata_sphinx_theme/layout.html#L118-L129
    "show_prev_next": False,
    "search_bar_text": "Search the docs ...",
    "navigation_with_keys": False,
    "collapse_navigation": False,
    "navigation_depth": 2,
    "show_nav_level": 1,
    "show_toc_level": 1,
    "navbar_align": "left",
    "header_links_before_dropdown": 5,
    "header_dropdown_text": "More",
    # The switcher requires a JSON file with the list of documentation versions, which
    # is generated by the script `build_tools/circle/list_versions.py` and placed under
    # the `js/` static directory; it will then be copied to the `_static` directory in
    # the built documentation
    "switcher": {
        "json_url": "https://scikit-learn.org/dev/_static/versions.json",
        "version_match": release,
    },
    # check_switcher may be set to False if docbuild pipeline fails. See
    # https://pydata-sphinx-theme.readthedocs.io/en/stable/user_guide/version-dropdown.html#configure-switcher-json-url
    "check_switcher": True,
    "pygments_light_style": "tango",
    "pygments_dark_style": "monokai",
    "logo": {
        "alt_text": "scikit-learn homepage",
        "image_relative": "logos/scikit-learn-logo-small.png",
        "image_light": "logos/scikit-learn-logo-small.png",
        "image_dark": "logos/scikit-learn-logo-small.png",
    },
    "surface_warnings": True,
    # -- Template placement in theme layouts ----------------------------------
    "navbar_start": ["navbar-logo"],
    # Note that the alignment of navbar_center is controlled by navbar_align
    "navbar_center": ["navbar-nav"],
    "navbar_end": ["theme-switcher", "navbar-icon-links", "version-switcher"],
    # navbar_persistent is persistent right (even when on mobiles)
    "navbar_persistent": ["search-button"],
    "article_header_start": ["breadcrumbs"],
    "article_header_end": [],
    "article_footer_items": ["prev-next"],
    "content_footer_items": [],
    # Use html_sidebars that map page patterns to list of sidebar templates
    "primary_sidebar_end": [],
    "footer_start": ["copyright"],
    "footer_center": [],
    "footer_end": [],
    # When specified as a dictionary, the keys should follow glob-style patterns, as in
    # https://www.sphinx-doc.org/en/master/usage/configuration.html#confval-exclude_patterns
    # In particular, "**" specifies the default for all pages
    # Use :html_theme.sidebar_secondary.remove: for file-wide removal
    "secondary_sidebar_items": {
        "**": [
            "page-toc",
            "sourcelink",
            # Sphinx-Gallery-specific sidebar components
            # https://sphinx-gallery.github.io/stable/advanced.html#using-sphinx-gallery-sidebar-components
            "sg_download_links",
            "sg_launcher_links",
        ],
    },
    "show_version_warning_banner": True,
    "announcement": (
        '<a href="https://forms.gle/zUXvWjGUN1nWhJ2V6">Help us make '
        "<code>scikit-learn</code> better! The 2024 user survey is now live.</a>"
    ),
}

# Add any paths that contain custom themes here, relative to this directory.
# html_theme_path = ["themes"]

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
# html_title = None

# A shorter title for the navigation bar.  Default is the same as html_title.
html_short_title = "scikit-learn"

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
html_favicon = "logos/favicon.ico"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["images", "css", "js"]

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
# html_last_updated_fmt = '%b %d, %Y'

# Custom sidebar templates, maps document names to template names.
# Workaround for removing the left sidebar on pages without TOC
# A better solution would be to follow the merge of:
# https://github.com/pydata/pydata-sphinx-theme/pull/1682
html_sidebars = {
    "install": [],
    "getting_started": [],
    "glossary": [],
    "faq": [],
    "support": [],
    "related_projects": [],
    "roadmap": [],
    "governance": [],
    "about": [],
}

# Additional templates that should be rendered to pages, maps page names to
# template names.
html_additional_pages = {"index": "index.html"}

# Additional files to copy
# html_extra_path = []

# Additional JS files
html_js_files = [
    "scripts/dropdown.js",
    "scripts/version-switcher.js",
]

# Compile scss files into css files using sphinxcontrib-sass
sass_src_dir, sass_out_dir = "scss", "css/styles"
sass_targets = {
    f"{file.stem}.scss": f"{file.stem}.css"
    for file in Path(sass_src_dir).glob("*.scss")
}

# Additional CSS files, should be subset of the values of `sass_targets`
html_css_files = ["styles/colors.css", "styles/custom.css"]


def add_js_css_files(app, pagename, templatename, context, doctree):
    """Load additional JS and CSS files only for certain pages.

    Note that `html_js_files` and `html_css_files` are included in all pages and
    should be used for the ones that are used by multiple pages. All page-specific
    JS and CSS files should be added here instead.
    """
    if pagename == "api/index":
        # External: jQuery and DataTables
        app.add_js_file("https://code.jquery.com/jquery-3.7.0.js")
        app.add_js_file("https://cdn.datatables.net/2.0.0/js/dataTables.min.js")
        app.add_css_file(
            "https://cdn.datatables.net/2.0.0/css/dataTables.dataTables.min.css"
        )
        # Internal: API search intialization and styling
        app.add_js_file("scripts/api-search.js")
        app.add_css_file("styles/api-search.css")
    elif pagename == "index":
        app.add_css_file("styles/index.css")
    elif pagename.startswith("modules/generated/"):
        app.add_css_file("styles/api.css")


# If false, no module index is generated.
html_domain_indices = False

# If false, no index is generated.
html_use_index = False

# If true, the index is split into individual pages for each letter.
# html_split_index = False

# If true, links to the reST sources are added to the pages.
# html_show_sourcelink = True

# If true, an OpenSearch description file will be output, and all pages will
# contain a <link> tag referring to it.  The value of this option must be the
# base URL from which the finished HTML is served.
# html_use_opensearch = ''

# If nonempty, this is the file name suffix for HTML files (e.g. ".xhtml").
# html_file_suffix = ''

# Output file base name for HTML help builder.
htmlhelp_basename = "scikit-learndoc"

# If true, the reST sources are included in the HTML build as _sources/name.
html_copy_source = True

# Adds variables into templates
html_context = {}
# finds latest release highlights and places it into HTML context for
# index.html
release_highlights_dir = Path("..") / "examples" / "release_highlights"
# Finds the highlight with the latest version number
latest_highlights = sorted(release_highlights_dir.glob("plot_release_highlights_*.py"))[
    -1
]
latest_highlights = latest_highlights.with_suffix("").name
html_context["release_highlights"] = (
    f"auto_examples/release_highlights/{latest_highlights}"
)

# get version from highlight name assuming highlights have the form
# plot_release_highlights_0_22_0
highlight_version = ".".join(latest_highlights.split("_")[-3:-1])
html_context["release_highlights_version"] = highlight_version


# redirects dictionary maps from old links to new links
redirects = {
    "documentation": "index",
    "contents": "index",
    "preface": "index",
    "modules/classes": "api/index",
    "tutorial/machine_learning_map/index": "machine_learning_map",
    "auto_examples/feature_selection/plot_permutation_test_for_classification": (
        "auto_examples/model_selection/plot_permutation_tests_for_classification"
    ),
    "modules/model_persistence": "model_persistence",
    "auto_examples/linear_model/plot_bayesian_ridge": (
        "auto_examples/linear_model/plot_ard"
    ),
    "auto_examples/model_selection/grid_search_text_feature_extraction": (
        "auto_examples/model_selection/plot_grid_search_text_feature_extraction"
    ),
    "auto_examples/datasets/plot_digits_last_image": (
        "auto_examples/exercises/plot_digits_classification_exercises"
    ),
    "auto_examples/miscellaneous/plot_changed_only_pprint_parameter": (
        "auto_examples/miscellaneous/plot_estimator_representation"
    ),
    "auto_examples/decomposition/plot_beta_divergence": (
        "auto_examples/applications/plot_topics_extraction_with_nmf_lda"
    ),
    "auto_examples/svm/plot_svm_nonlinear": "auto_examples/svm/plot_svm_kernels",
    "auto_examples/ensemble/plot_adaboost_hastie_10_2": (
        "auto_examples/ensemble/plot_adaboost_multiclass"
    ),
    "auto_examples/decomposition/plot_pca_3d": (
        "auto_examples/decomposition/plot_pca_iris"
    ),
    "auto_examples/exercises/plot_cv_digits": (
        "auto_examples/model_selection/plot_nested_cross_validation_iris"
    ),
    "auto_examples/linear_model/plot_lasso_lars": (
        "auto_examples/linear_model/plot_lasso_lasso_lars_elasticnet_path"
    ),
    "auto_examples/linear_model/plot_lasso_coordinate_descent_path": (
        "auto_examples/linear_model/plot_lasso_lasso_lars_elasticnet_path"
    ),
}
html_context["redirects"] = redirects
for old_link in redirects:
    html_additional_pages[old_link] = "redirects.html"

# See https://github.com/scikit-learn/scikit-learn/pull/22550
html_context["is_devrelease"] = parsed_version.is_devrelease


# -- Options for LaTeX output ------------------------------------------------
latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    # 'papersize': 'letterpaper',
    # The font size ('10pt', '11pt' or '12pt').
    # 'pointsize': '10pt',
    # Additional stuff for the LaTeX preamble.
    "preamble": r"""
        \usepackage{amsmath}\usepackage{amsfonts}\usepackage{bm}
        \usepackage{morefloats}\usepackage{enumitem} \setlistdepth{10}
        \let\oldhref\href
        \renewcommand{\href}[2]{\oldhref{#1}{\hbox{#2}}}
        """
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, documentclass
# [howto/manual]).
latex_documents = [
    (
        "contents",
        "user_guide.tex",
        "scikit-learn user guide",
        "scikit-learn developers",
        "manual",
    ),
]

# The name of an image file (relative to this directory) to place at the top of
# the title page.
latex_logo = "logos/scikit-learn-logo.png"

# Documents to append as an appendix to all manuals.
# latex_appendices = []

# If false, no module index is generated.
latex_domain_indices = False

trim_doctests_flags = True

# intersphinx configuration
intersphinx_mapping = {
    "python": ("https://docs.python.org/{.major}".format(sys.version_info), None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "matplotlib": ("https://matplotlib.org/", None),
    "pandas": ("https://pandas.pydata.org/pandas-docs/stable/", None),
    "joblib": ("https://joblib.readthedocs.io/en/latest/", None),
    "seaborn": ("https://seaborn.pydata.org/", None),
    "skops": ("https://skops.readthedocs.io/en/stable/", None),
}

v = parse(release)
if v.release is None:
    raise ValueError(
        "Ill-formed version: {!r}. Version should follow PEP440".format(version)
    )

if v.is_devrelease:
    binder_branch = "main"
else:
    major, minor = v.release[:2]
    binder_branch = "{}.{}.X".format(major, minor)


class SubSectionTitleOrder:
    """Sort example gallery by title of subsection.

    Assumes README.txt exists for all subsections and uses the subsection with
    dashes, '---', as the adornment.
    """

    def __init__(self, src_dir):
        self.src_dir = src_dir
        self.regex = re.compile(r"^([\w ]+)\n-", re.MULTILINE)

    def __repr__(self):
        return "<%s>" % (self.__class__.__name__,)

    def __call__(self, directory):
        src_path = os.path.normpath(os.path.join(self.src_dir, directory))

        # Forces Release Highlights to the top
        if os.path.basename(src_path) == "release_highlights":
            return "0"

        readme = os.path.join(src_path, "README.txt")

        try:
            with open(readme, "r") as f:
                content = f.read()
        except FileNotFoundError:
            return directory

        title_match = self.regex.search(content)
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
        code_lines.append("%pip install plotly")
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

    if code_lines:
        code_lines = ["# JupyterLite-specific code"] + code_lines
        code = "\n".join(code_lines)
        add_code_cell(dummy_notebook_content, code)

    notebook_content["cells"] = (
        dummy_notebook_content["cells"] + notebook_content["cells"]
    )


default_global_config = sklearn.get_config()


def reset_sklearn_config(gallery_conf, fname):
    """Reset sklearn config to default values."""
    sklearn.set_config(**default_global_config)


sg_examples_dir = "../examples"
sg_gallery_dir = "auto_examples"
sphinx_gallery_conf = {
    "doc_module": "sklearn",
    "backreferences_dir": os.path.join("modules", "generated"),
    "show_memory": False,
    "reference_url": {"sklearn": None},
    "examples_dirs": [sg_examples_dir],
    "gallery_dirs": [sg_gallery_dir],
    "subsection_order": SubSectionTitleOrder(sg_examples_dir),
    "within_subsection_order": SKExampleTitleSortKey,
    "binder": {
        "org": "scikit-learn",
        "repo": "scikit-learn",
        "binderhub_url": "https://mybinder.org",
        "branch": binder_branch,
        "dependencies": "./binder/requirements.txt",
        "use_jupyter_lab": True,
    },
    # avoid generating too many cross links
    "inspect_global_variables": False,
    "remove_config_comments": True,
    "plot_gallery": "True",
    "recommender": {"enable": True, "n_examples": 4, "min_df": 12},
    "reset_modules": ("matplotlib", "seaborn", reset_sklearn_config),
}
if with_jupyterlite:
    sphinx_gallery_conf["jupyterlite"] = {
        "notebook_modification_function": notebook_modification_function
    }

# For the index page of the gallery and each nested section, we hide the secondary
# sidebar by specifying an empty list (no components), because there is no meaningful
# in-page toc for these pages, and they are generated so "sourcelink" is not useful
# either.
html_theme_options["secondary_sidebar_items"][f"{sg_gallery_dir}/index"] = []
for sub_sg_dir in (Path(".") / sg_examples_dir).iterdir():
    if sub_sg_dir.is_dir():
        html_theme_options["secondary_sidebar_items"][
            f"{sg_gallery_dir}/{sub_sg_dir.name}/index"
        ] = []


# The following dictionary contains the information used to create the
# thumbnails for the front page of the scikit-learn home page.
# key: first image in set
# values: (number of plot in set, height of thumbnail)
carousel_thumbs = {"sphx_glr_plot_classifier_comparison_001.png": 600}


# enable experimental module so that experimental estimators can be
# discovered properly by sphinx
from sklearn.experimental import enable_iterative_imputer  # noqa
from sklearn.experimental import enable_halving_search_cv  # noqa


def make_carousel_thumbs(app, exception):
    """produces the final resized carousel images"""
    if exception is not None:
        return
    print("Preparing carousel images")

    image_dir = os.path.join(app.builder.outdir, "_images")
    for glr_plot, max_width in carousel_thumbs.items():
        image = os.path.join(image_dir, glr_plot)
        if os.path.exists(image):
            c_thumb = os.path.join(image_dir, glr_plot[:-4] + "_carousel.png")
            sphinx_gallery.gen_rst.scale_image(image, c_thumb, max_width, 190)


def filter_search_index(app, exception):
    if exception is not None:
        return

    # searchindex only exist when generating html
    if app.builder.name != "html":
        return

    print("Removing methods from search index")

    searchindex_path = os.path.join(app.builder.outdir, "searchindex.js")
    with open(searchindex_path, "r") as f:
        searchindex_text = f.read()

    searchindex_text = re.sub(r"{__init__.+?}", "{}", searchindex_text)
    searchindex_text = re.sub(r"{__call__.+?}", "{}", searchindex_text)

    with open(searchindex_path, "w") as f:
        f.write(searchindex_text)


# Config for sphinx_issues

# we use the issues path for PRs since the issues URL will forward
issues_github_path = "scikit-learn/scikit-learn"


def disable_plot_gallery_for_linkcheck(app):
    if app.builder.name == "linkcheck":
        sphinx_gallery_conf["plot_gallery"] = "False"


def setup(app):
    # do not run the examples when using linkcheck by using a small priority
    # (default priority is 500 and sphinx-gallery using builder-inited event too)
    app.connect("builder-inited", disable_plot_gallery_for_linkcheck, priority=50)

    # triggered just before the HTML for an individual page is created
    app.connect("html-page-context", add_js_css_files)

    # to hide/show the prompt in code examples
    app.connect("build-finished", make_carousel_thumbs)
    app.connect("build-finished", filter_search_index)


# The following is used by sphinx.ext.linkcode to provide links to github
linkcode_resolve = make_linkcode_resolve(
    "sklearn",
    (
        "https://github.com/scikit-learn/"
        "scikit-learn/blob/{revision}/"
        "{package}/{path}#L{lineno}"
    ),
)

warnings.filterwarnings(
    "ignore",
    category=UserWarning,
    message=(
        "Matplotlib is currently using agg, which is a"
        " non-GUI backend, so cannot show the figure."
    ),
)
if os.environ.get("SKLEARN_WARNINGS_AS_ERRORS", "0") != "0":
    turn_warnings_into_errors()

# maps functions with a class name that is indistinguishable when case is
# ignore to another filename
autosummary_filename_map = {
    "sklearn.cluster.dbscan": "dbscan-function",
    "sklearn.covariance.oas": "oas-function",
    "sklearn.decomposition.fastica": "fastica-function",
}


# Config for sphinxext.opengraph

ogp_site_url = "https://scikit-learn/stable/"
ogp_image = "https://scikit-learn.org/stable/_static/scikit-learn-logo-small.png"
ogp_use_first_image = True
ogp_site_name = "scikit-learn"

# Config for linkcheck that checks the documentation for broken links

# ignore all links in 'whats_new' to avoid doing many github requests and
# hitting the github rate threshold that makes linkcheck take a lot of time
linkcheck_exclude_documents = [r"whats_new/.*"]

# default timeout to make some sites links fail faster
linkcheck_timeout = 10

# Allow redirects from doi.org
linkcheck_allowed_redirects = {r"https://doi.org/.+": r".*"}
linkcheck_ignore = [
    # ignore links to local html files e.g. in image directive :target: field
    r"^..?/",
    # ignore links to specific pdf pages because linkcheck does not handle them
    # ('utf-8' codec can't decode byte error)
    r"http://www.utstat.toronto.edu/~rsalakhu/sta4273/notes/Lecture2.pdf#page=.*",
    (
        "https://www.fordfoundation.org/media/2976/roads-and-bridges"
        "-the-unseen-labor-behind-our-digital-infrastructure.pdf#page=.*"
    ),
    # links falsely flagged as broken
    (
        "https://www.researchgate.net/publication/"
        "233096619_A_Dendrite_Method_for_Cluster_Analysis"
    ),
    (
        "https://www.researchgate.net/publication/221114584_Random_Fourier"
        "_Approximations_for_Skewed_Multiplicative_Histogram_Kernels"
    ),
    (
        "https://www.researchgate.net/publication/4974606_"
        "Hedonic_housing_prices_and_the_demand_for_clean_air"
    ),
    (
        "https://www.researchgate.net/profile/Anh-Huy-Phan/publication/220241471_Fast_"
        "Local_Algorithms_for_Large_Scale_Nonnegative_Matrix_and_Tensor_Factorizations"
    ),
    "https://doi.org/10.13140/RG.2.2.35280.02565",
    (
        "https://www.microsoft.com/en-us/research/uploads/prod/2006/01/"
        "Bishop-Pattern-Recognition-and-Machine-Learning-2006.pdf"
    ),
    "https://www.microsoft.com/en-us/research/wp-content/uploads/2016/02/tr-99-87.pdf",
    "https://microsoft.com/",
    "https://www.jstor.org/stable/2984099",
    "https://stat.uw.edu/sites/default/files/files/reports/2000/tr371.pdf",
    # Broken links from testimonials
    "http://www.bestofmedia.com",
    "http://www.data-publica.com/",
    "https://livelovely.com",
    "https://www.mars.com/global",
    "https://www.yhat.com",
    # Ignore some dynamically created anchors. See
    # https://github.com/sphinx-doc/sphinx/issues/9016 for more details about
    # the github example
    r"https://github.com/conda-forge/miniforge#miniforge",
    r"https://github.com/joblib/threadpoolctl/"
    "#setting-the-maximum-size-of-thread-pools",
    r"https://stackoverflow.com/questions/5836335/"
    "consistently-create-same-random-numpy-array/5837352#comment6712034_5837352",
]

# Config for sphinx-remove-toctrees

remove_from_toctrees = ["metadata_routing.rst"]

# Use a browser-like user agent to avoid some "403 Client Error: Forbidden for
# url" errors. This is taken from the variable navigator.userAgent inside a
# browser console.
user_agent = (
    "Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:100.0) Gecko/20100101 Firefox/100.0"
)

# Use Github token from environment variable to avoid Github rate limits when
# checking Github links
github_token = os.getenv("GITHUB_TOKEN")

if github_token is None:
    linkcheck_request_headers = {}
else:
    linkcheck_request_headers = {
        "https://github.com/": {"Authorization": f"token {github_token}"},
    }


def infer_next_release_versions():
    """Infer the most likely next release versions to make."""
    all_version_full = {"rc": "0.99.0rc1", "final": "0.99.0", "bf": "0.98.1"}
    all_version_short = {"rc": "0.99", "final": "0.99", "bf": "0.98"}
    all_previous_tag = {"rc": "unused", "final": "0.98.33", "bf": "0.97.22"}

    try:
        # Fetch the version switcher JSON; see `html_theme_options` for more details
        versions_json = json.loads(
            urlopen(html_theme_options["switcher"]["json_url"], timeout=10).read()
        )

        # See `build_tools/circle/list_versions.py`, stable is always the second entry
        stable_version = parse(versions_json[1]["version"])
        last_stable_version = parse(versions_json[2]["version"])
        next_major_minor = f"{stable_version.major}.{stable_version.minor + 1}"

        # RC
        all_version_full["rc"] = f"{next_major_minor}.0rc1"
        all_version_short["rc"] = next_major_minor

        # Major/Minor final
        all_version_full["final"] = f"{next_major_minor}.0"
        all_version_short["final"] = next_major_minor
        all_previous_tag["final"] = stable_version.base_version

        # Bug-fix
        all_version_full["bf"] = (
            f"{stable_version.major}.{stable_version.minor}.{stable_version.micro + 1}"
        )
        all_version_short["bf"] = f"{stable_version.major}.{stable_version.minor}"
        all_previous_tag["bf"] = last_stable_version.base_version
    except Exception as e:
        logger.warning(
            "Failed to infer all possible next release versions because of "
            f"{type(e).__name__}: {e}"
        )

    return {
        "version_full": all_version_full,
        "version_short": all_version_short,
        "previous_tag": all_previous_tag,
    }


# -- Convert .rst.template files to .rst ---------------------------------------

from api_reference import API_REFERENCE, DEPRECATED_API_REFERENCE

from sklearn._min_dependencies import dependent_packages

# If development build, link to local page in the top navbar; otherwise link to the
# development version; see https://github.com/scikit-learn/scikit-learn/pull/22550
if parsed_version.is_devrelease:
    development_link = "developers/index"
else:
    development_link = "https://scikit-learn.org/dev/developers/index.html"

# Define the templates and target files for conversion
# Each entry is in the format (template name, file name, kwargs for rendering)
rst_templates = [
    ("index", "index", {"development_link": development_link}),
    (
        "developers/maintainer",
        "developers/maintainer",
        {"inferred": infer_next_release_versions()},
    ),
    (
        "min_dependency_table",
        "min_dependency_table",
        {"dependent_packages": dependent_packages},
    ),
    (
        "min_dependency_substitutions",
        "min_dependency_substitutions",
        {"dependent_packages": dependent_packages},
    ),
    (
        "api/index",
        "api/index",
        {
            "API_REFERENCE": sorted(API_REFERENCE.items(), key=lambda x: x[0]),
            "DEPRECATED_API_REFERENCE": sorted(
                DEPRECATED_API_REFERENCE.items(), key=lambda x: x[0], reverse=True
            ),
        },
    ),
]

# Convert each module API reference page
for module in API_REFERENCE:
    rst_templates.append(
        (
            "api/module",
            f"api/{module}",
            {"module": module, "module_info": API_REFERENCE[module]},
        )
    )

# Convert the deprecated API reference page (if there exists any)
if DEPRECATED_API_REFERENCE:
    rst_templates.append(
        (
            "api/deprecated",
            "api/deprecated",
            {
                "DEPRECATED_API_REFERENCE": sorted(
                    DEPRECATED_API_REFERENCE.items(), key=lambda x: x[0], reverse=True
                )
            },
        )
    )

for rst_template_name, rst_target_name, kwargs in rst_templates:
    # Read the corresponding template file into jinja2
    with (Path(".") / f"{rst_template_name}.rst.template").open(
        "r", encoding="utf-8"
    ) as f:
        t = jinja2.Template(f.read())

    # Render the template and write to the target
    with (Path(".") / f"{rst_target_name}.rst").open("w", encoding="utf-8") as f:
        f.write(t.render(**kwargs))
