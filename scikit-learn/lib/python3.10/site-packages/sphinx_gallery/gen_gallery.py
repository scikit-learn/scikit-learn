# Author: Óscar Nájera
# License: 3-clause BSD
"""Sphinx-Gallery Generator.

Attaches Sphinx-Gallery to Sphinx in order to generate the galleries
when building the documentation.
"""

import codecs
import copy
import os
import pathlib
import re
from datetime import datetime, timedelta
from difflib import get_close_matches
from itertools import chain
from pathlib import Path
from textwrap import indent
from xml.sax.saxutils import escape, quoteattr

import sphinx.util
from docutils import nodes
from sphinx.errors import ConfigError, ExtensionError
from sphinx.util.console import blue, bold, purple, red

from . import __version__ as _sg_version
from . import glr_path_static
from .backreferences import _finalize_backreferences
from .directives import ImageSg, MiniGallery, imagesg_addnode
from .docs_resolv import embed_code_links
from .downloads import generate_zipfiles
from .gen_rst import (
    SPHX_GLR_SIG,
    _get_call_memory_and_base,
    _get_callables,
    _get_gallery_header,
    generate_dir_rst,
)
from .interactive_example import (
    check_binder_conf,
    check_jupyterlite_conf,
    copy_binder_files,
    create_jupyterlite_contents,
    post_configure_jupyterlite_sphinx,
    pre_configure_jupyterlite_sphinx,
)
from .recommender import ExampleRecommender, _write_recommendations
from .scrapers import _import_matplotlib
from .sorting import ExplicitOrder
from .utils import (
    _collect_gallery_files,
    _combine_backreferences,
    _format_toctree,
    _has_graphviz,
    _has_optipng,
    _has_pypandoc,
    _read_json,
    _replace_md5,
    _write_json,
)

_KNOWN_CSS = (
    "sg_gallery",
    "sg_gallery-binder",
    "sg_gallery-dataframe",
    "sg_gallery-rendered-html",
)


class DefaultResetArgv:
    """Provides default 'reset_argv' callable that returns empty list."""

    def __repr__(self):
        return "DefaultResetArgv"

    def __call__(self, gallery_conf, script_vars):
        """Return empty list."""
        return []


DEFAULT_GALLERY_CONF = {
    "filename_pattern": re.escape(os.sep) + "plot",
    "ignore_pattern": r"__init__\.py",
    "examples_dirs": os.path.join("..", "examples"),
    "example_extensions": {".py"},
    "filetype_parsers": {},
    "notebook_extensions": {".py"},
    "reset_argv": DefaultResetArgv(),
    "subsection_order": None,
    "within_subsection_order": "NumberOfCodeLinesSortKey",
    "minigallery_sort_order": None,
    "gallery_dirs": "auto_examples",
    "backreferences_dir": None,
    "doc_module": (),
    "exclude_implicit_doc": set(),
    "reference_url": {},
    "capture_repr": ("_repr_html_", "__repr__"),
    "ignore_repr_types": r"",
    # Build options
    # -------------
    # 'plot_gallery' also accepts strings that evaluate to a bool, e.g. "True",
    # "False", "1", "0" so that they can be easily set via command line
    # switches of sphinx-build
    "plot_gallery": "True",
    "download_all_examples": True,
    "abort_on_example_error": False,
    "only_warn_on_example_error": False,
    "recommender": {"enable": False},
    "failing_examples": {},
    "passing_examples": [],
    "stale_examples": [],  # ones that did not need to be run due to md5sum
    "run_stale_examples": False,
    "expected_failing_examples": set(),
    "thumbnail_size": (400, 280),  # Default CSS does 0.4 scaling (160, 112)
    "min_reported_time": 0,
    "write_computation_times": os.getenv("SOURCE_DATE_EPOCH") is None,
    "binder": {},
    "jupyterlite": {},
    "promote_jupyter_magic": False,
    "image_scrapers": ("matplotlib",),
    "compress_images": (),
    "reset_modules": ("matplotlib", "seaborn"),
    "reset_modules_order": "before",
    "first_notebook_cell": None,
    "last_notebook_cell": None,
    "notebook_images": False,
    "pypandoc": False,
    "remove_config_comments": False,
    "show_memory": False,
    "show_signature": True,
    "junit": "",
    "log_level": {"backreference_missing": "warning"},
    "inspect_global_variables": True,
    "css": _KNOWN_CSS,
    "matplotlib_animations": False,
    "image_srcset": [],
    "default_thumb_file": None,
    "line_numbers": False,
    "nested_sections": True,
    "prefer_full_module": set(),
    "api_usage_ignore": ".*__.*__",
    "show_api_usage": False,  # if this changes, change write_api_entries, too
    "copyfile_regex": "",
    "parallel": False,
}

logger = sphinx.util.logging.getLogger("sphinx-gallery")


def _bool_eval(x):
    """Evaluate bool only configs, to allow setting via -D on the command line."""
    if isinstance(x, str):
        try:
            x = eval(x)
        except TypeError:
            pass
    return bool(x)


def _update_gallery_conf_exclude_implicit_doc(gallery_conf):
    """Update gallery config exclude_implicit_doc.

    This is separate function for better testability.
    """
    # prepare regex for exclusions from implicit documentation, ensuring that what
    # gets complied has a stable __repr__ (i.e., by sorting the exclude_implicit_doc
    # set before joining)
    exclude_regex = (
        re.compile("|".join(sorted(gallery_conf["exclude_implicit_doc"])))
        if gallery_conf["exclude_implicit_doc"]
        else False
    )
    gallery_conf["exclude_implicit_doc_regex"] = exclude_regex


def _update_gallery_conf_builder_inited(
    sphinx_gallery_conf,
    src_dir,
    plot_gallery=True,
    abort_on_example_error=False,
    builder_name="html",
):
    sphinx_gallery_conf.update(plot_gallery=plot_gallery)
    sphinx_gallery_conf.update(abort_on_example_error=abort_on_example_error)
    sphinx_gallery_conf["src_dir"] = src_dir
    # Make it easy to know which builder we're in
    sphinx_gallery_conf["builder_name"] = builder_name


def _check_extra_config_keys(gallery_conf, sphinx_gallery_conf, check_keys):
    """Check SG config keys, optionally raising if any extra keys present."""
    options = sorted(gallery_conf)
    extra_keys = sorted(set(sphinx_gallery_conf) - set(options))
    if extra_keys and check_keys:
        msg = "Unknown key(s) in sphinx_gallery_conf:\n"
        for key in extra_keys:
            options = get_close_matches(key, options, cutoff=0.66)
            msg += repr(key)
            if len(options) == 1:
                msg += f", did you mean {options[0]!r}?"
            elif len(options) > 1:
                msg += f", did you mean one of {options!r}?"
            msg += "\n"
        raise ConfigError(msg.strip())


def _check_config_type(
    gallery_conf,
    conf_key,
    types,
    str_to_list=False,
    allow_none=False,
):
    """Check config type, optionally converting str to list or allowing None."""
    conf_value = gallery_conf[conf_key]
    # Early exit if type correct
    if isinstance(conf_value, types) or (allow_none and conf_value is None):
        if str_to_list:
            if isinstance(conf_value, str):
                gallery_conf[conf_key] = [conf_value]
            return
        return
    # Otherwise raise error
    msg = "'{conf_key}' config allowed types: {types}{or_none}. Got {conf_type}."
    if isinstance(types, type):
        types = (types,)
    str_types = [t.__name__ for t in types]
    or_none = " or None" if allow_none else ""
    conf_type = type(conf_value).__name__
    msg = msg.format(
        conf_key=conf_key,
        types=str_types,
        or_none=or_none,
        conf_type=conf_type,
    )
    raise ConfigError(msg)


def _check_image_srcset(gallery_conf):
    """Check `_check_image_srcset`, convert to float and removing '1'."""
    _check_config_type(gallery_conf, "image_srcset", (list, tuple), str_to_list=True)
    srcset_mult_facs = set()
    for st in gallery_conf["image_srcset"]:
        if not (isinstance(st, str) and st[-1:] == "x"):
            raise ConfigError(
                f"Invalid value for image_srcset parameter: {st!r}. "
                "Must be a list of strings with the multiplicative "
                'factor followed by an "x".  e.g. ["2.0x", "1.5x"]'
            )
        # "2x" -> "2.0"
        srcset_mult_facs.add(float(st[:-1]))
    srcset_mult_facs -= {1}  # 1x is always saved.
    gallery_conf["image_srcset"] = [*sorted(srcset_mult_facs)]


def _check_compress_images(gallery_conf):
    """Check `compress_images`, getting any command line args."""
    _check_config_type(
        gallery_conf,
        "compress_images",
        (str, tuple, list),
        str_to_list=True,
    )
    compress_images = gallery_conf["compress_images"]
    compress_images = list(compress_images)
    allowed_values = ("images", "thumbnails")
    pops = list()
    # Get command-line switches
    for ki, kind in enumerate(compress_images):
        if kind not in allowed_values:
            if kind.startswith("-"):
                pops.append(ki)
                continue
            raise ConfigError(
                "All entries in compress_images must be one of "
                f"{allowed_values} or a command-line switch "
                f'starting with "-", got {kind!r}'
            )
    compress_images_args = [compress_images.pop(p) for p in pops[::-1]]
    if len(compress_images) and not _has_optipng():
        logger.warning(
            "optipng binaries not found, PNG %s will not be optimized",
            " and ".join(compress_images),
        )
        compress_images = ()
    gallery_conf["compress_images"] = compress_images
    gallery_conf["compress_images_args"] = compress_images_args


def _check_matplotlib_animations(gallery_conf, app):
    """Check `matplotlib_animations` config."""
    animations = gallery_conf["matplotlib_animations"]
    if isinstance(animations, bool):
        # We allow single boolean for backwards compatibility reasons
        animations = (animations,)
        fmt = None
    if not (len(animations) in (1, 2) and isinstance(enabled := animations[0], bool)):
        raise ConfigError(
            "'matplotlib_animations' must be a single bool or "
            f"(enabled: bool, format: str), not {animations!r}"
        )
    # Handle file format
    if len(animations) > 1:
        fmt = animations[1]
        if fmt is not None:
            if not isinstance(fmt, str):
                raise ConfigError(
                    "'matplotlib_animations' file format must be a string or None"
                )
            if fmt not in ("html5", "jshtml"):
                if app is not None:
                    # Other formats mean animations saved externally and require
                    # this `video` extension to embed them into the HTML
                    try:
                        app.setup_extension("sphinxcontrib.video")
                    except ExtensionError as e:
                        raise ConfigError(
                            f"'matplotlib_animations' specifies file format: {fmt}; "
                            f"this requires the sphinxcontrib.video package."
                        ) from e

    gallery_conf["matplotlib_animations"] = (enabled, fmt)


def _check_pypandoc_config(gallery_conf):
    """Check `pypandoc` config."""
    pypandoc = gallery_conf["pypandoc"]
    _check_config_type(gallery_conf, "pypandoc", (dict, bool))

    gallery_conf["pypandoc"] = dict() if pypandoc is True else pypandoc
    has_pypandoc, version = _has_pypandoc()
    if isinstance(gallery_conf["pypandoc"], dict) and has_pypandoc is None:
        logger.warning(
            "'pypandoc' not available. Using Sphinx-Gallery to "
            "convert rst text blocks to markdown for .ipynb files."
        )
        gallery_conf["pypandoc"] = False
    elif isinstance(gallery_conf["pypandoc"], dict):
        logger.info(
            "Using pandoc version: %s to convert rst text blocks to "
            "markdown for .ipynb files",
            version,
        )
    else:
        logger.info(
            "Using Sphinx-Gallery to convert rst text blocks to "
            "markdown for .ipynb files."
        )
    if isinstance(pypandoc, dict):
        accepted_keys = ("extra_args", "filters")
        for key in pypandoc:
            if key not in accepted_keys:
                raise ConfigError(
                    "'pypandoc' only accepts the following key "
                    f"values: {accepted_keys}, got: {key}."
                )


def _fill_gallery_conf_defaults(sphinx_gallery_conf, app=None, check_keys=True):
    """Handle user configs, update default gallery configs and check values."""
    gallery_conf = copy.deepcopy(DEFAULT_GALLERY_CONF)
    _check_extra_config_keys(gallery_conf, sphinx_gallery_conf, check_keys)

    gallery_conf.update(sphinx_gallery_conf)
    # XXX anything that can only be a bool (rather than str) should probably be
    # evaluated this way as it allows setting via -D on the command line
    for key in (
        "promote_jupyter_magic",
        "run_stale_examples",
    ):
        gallery_conf[key] = _bool_eval(gallery_conf[key])

    gallery_conf["default_role"] = ""
    gallery_conf["source_suffix"] = {".rst": "restructuredtext"}
    if app is not None:
        if app.config["default_role"]:
            gallery_conf["default_role"] = app.config["default_role"]
        gallery_conf["source_suffix"] = app.config["source_suffix"]
        if isinstance(gallery_conf["source_suffix"], str):
            gallery_conf["source_suffix"] = {gallery_conf["source_suffix"]: None}

    # Check capture_repr
    _check_config_type(gallery_conf, "capture_repr", (tuple, list), str_to_list=True)
    supported_reprs = ["__repr__", "__str__", "_repr_html_"]
    for rep in gallery_conf["capture_repr"]:
        if rep not in supported_reprs:
            raise ConfigError(
                "All entries in 'capture_repr' must be one "
                f"of {supported_reprs}, got: {rep}"
            )
    # Check ignore_repr_types
    _check_config_type(gallery_conf, "ignore_repr_types", str)

    # Check parallel
    _check_config_type(gallery_conf, "parallel", (bool, int))
    if gallery_conf["parallel"] is True:
        gallery_conf["parallel"] = app.parallel
    if gallery_conf["parallel"] == 1:
        gallery_conf["parallel"] = False
    if gallery_conf["parallel"]:
        try:
            import joblib  # noqa
        except Exception:
            raise ValueError("joblib must be importable when parallel mode is enabled")

    # deal with show_memory
    _get_call_memory_and_base(gallery_conf, update=True)

    # check callables
    for key in (
        "image_scrapers",
        "reset_argv",
        "minigallery_sort_order",
        "reset_modules",
    ):
        if key == "minigallery_sort_order" and gallery_conf[key] is None:
            continue
        _get_callables(gallery_conf, key)

    # Here we try to set up matplotlib but don't raise an error,
    # we will raise an error later when we actually try to use it
    # (if we do so) in scrapers.py.
    # In principle we could look to see if there is a matplotlib scraper
    # in our scrapers list, but this would be backward incompatible with
    # anyone using or relying on our Agg-setting behavior (e.g., for some
    # custom matplotlib SVG scraper as in our docs).
    # Eventually we can make this a config var like matplotlib_agg or something
    # if people need us not to set it to Agg.
    try:
        _import_matplotlib()
    except (ImportError, ValueError):
        pass

    # Check for srcset hidpi images
    _check_image_srcset(gallery_conf)
    # Check `compress_images`
    _check_compress_images(gallery_conf)
    # Check `matplotlib_animations`
    _check_matplotlib_animations(gallery_conf, app)

    # Check resetters
    _get_callables(gallery_conf, "reset_modules")

    _check_config_type(gallery_conf, "reset_modules_order", str)
    if gallery_conf["reset_modules_order"] not in ["before", "after", "both"]:
        raise ConfigError(
            "reset_modules_order must be in"
            "['before', 'after', 'both'], "
            f"got {gallery_conf['reset_modules_order']!r}"
        )

    # Ensure the first/last cell text is a string if we have it
    cell_config_keys = ("first_notebook_cell", "last_notebook_cell")
    for conf_key in cell_config_keys:
        _check_config_type(gallery_conf, conf_key, str, allow_none=True)

    # Check pypandoc
    _check_pypandoc_config(gallery_conf)

    gallery_conf["titles"] = {}
    # Ensure 'backreferences_dir' is str, pathlib.Path or None
    _check_config_type(
        gallery_conf,
        "backreferences_dir",
        (str, pathlib.Path),
        allow_none=True,
    )
    # if 'backreferences_dir' is pathlib.Path, make str for Python <=3.5
    # compatibility
    backref = gallery_conf["backreferences_dir"]
    if isinstance(backref, pathlib.Path):
        gallery_conf["backreferences_dir"] = str(backref)

    # binder
    gallery_conf["binder"] = check_binder_conf(gallery_conf["binder"])

    # jupyterlite
    gallery_conf["jupyterlite"] = check_jupyterlite_conf(
        gallery_conf["jupyterlite"],
        app,
    )

    # css
    _check_config_type(gallery_conf, "css", (list, tuple), str_to_list=True)
    for css in gallery_conf["css"]:
        if css not in _KNOWN_CSS:
            raise ConfigError(f"Unknown css {css!r}, must be one of {_KNOWN_CSS!r}")
        if app is not None:  # can be None in testing
            app.add_css_file(css + ".css")

    # check API usage
    _check_config_type(gallery_conf, "api_usage_ignore", str)
    if (
        not isinstance(gallery_conf["show_api_usage"], bool)
        and gallery_conf["show_api_usage"] != "unused"
    ):
        raise ConfigError(
            'gallery_conf["show_api_usage"] must be True, False or "unused", '
            f"got {gallery_conf['show_api_usage']}"
        )

    # check `within_subsection_order`
    _get_callables(
        gallery_conf, "within_subsection_order", src_dir=""
    )  # make sure it works

    _update_gallery_conf_exclude_implicit_doc(gallery_conf)

    return gallery_conf


def get_subsections(srcdir, examples_dir, gallery_conf, check_for_header=True):
    """Return the list of subsections of a gallery.

    Parameters
    ----------
    srcdir : str
        absolute path to directory containing conf.py
    examples_dir : str
        path to the examples directory relative to conf.py
    gallery_conf : Dict[str, Any]
        Sphinx-Gallery configuration dictionary.
    check_for_header : bool
        only return subfolders that contain a GALLERY_HEADER file, default True

    Returns
    -------
    out : list
        sorted list of gallery subsection folder names
    """
    if gallery_conf["subsection_order"] is None:
        sortkey = None
    else:
        subsec_order = gallery_conf.get("subsection_order")
        if isinstance(subsec_order, list):
            gallery_conf["subsection_order"] = ExplicitOrder(subsec_order)
        (sortkey,) = _get_callables(gallery_conf, "subsection_order")
    subfolders = [subfolder for subfolder in os.listdir(examples_dir)]
    if check_for_header:
        subfolders = [
            subfolder
            for subfolder in subfolders
            # Return is not `None` only when a gallery head file is found
            if _get_gallery_header(
                os.path.join(examples_dir, subfolder), gallery_conf, raise_error=False
            )
            is not None
        ]
    else:
        # just make sure its a directory, that is not `__pycache__`
        subfolders = [
            subfolder
            for subfolder in subfolders
            if (
                subfolder != "__pycache__"
                and os.path.isdir(os.path.join(examples_dir, subfolder))
            )
        ]

    base_examples_dir_path = os.path.relpath(examples_dir, srcdir)
    subfolders_with_path = [
        os.path.join(base_examples_dir_path, item) for item in subfolders
    ]
    sorted_subfolders = sorted(subfolders_with_path, key=sortkey)

    return [
        subfolders[i]
        for i in [subfolders_with_path.index(item) for item in sorted_subfolders]
    ]


def _prepare_sphx_glr_dirs(gallery_conf, srcdir):
    """Creates necessary folders for sphinx_gallery files."""
    examples_dirs = gallery_conf["examples_dirs"]
    gallery_dirs = gallery_conf["gallery_dirs"]

    if not isinstance(examples_dirs, list):
        examples_dirs = [examples_dirs]

    if not isinstance(gallery_dirs, list):
        gallery_dirs = [gallery_dirs]

    if len(examples_dirs) != len(gallery_dirs):
        logger.warning(
            "'examples_dirs' and 'gallery_dirs' are of different lengths. "
            "Surplus entries will be ignored."
        )

    if bool(gallery_conf["backreferences_dir"]):
        backreferences_dir = os.path.join(srcdir, gallery_conf["backreferences_dir"])
        os.makedirs(backreferences_dir, exist_ok=True)

    return list(zip(examples_dirs, gallery_dirs))


def _filter_tags(subsection_index_content):
    """Filter out tags from `subsection_index_content`."""
    tag_regex = r"^\.\.(\s+)\_(.+)\:(\s*)$"
    subsection_index_content = "\n".join(
        [
            line
            for line in subsection_index_content.splitlines()
            if re.match(tag_regex, line) is None
        ]
        + [""]
    )
    return subsection_index_content


def _finish_index_rst(
    app,
    gallery_conf,
    indexst,
    sg_root_index,
    subsection_index_files,
    gallery_dir_abs_path,
):
    """Add toctree, download and signature, if req, to index and write file."""
    # Generate toctree containing subsection index files
    if (
        sg_root_index
        and gallery_conf["nested_sections"] is True
        and len(subsection_index_files) > 0
    ):
        subsections_toctree = _format_toctree(
            subsection_index_files, includehidden=True
        )
        indexst += subsections_toctree

    # Always generate download zipfiles, only add to index.rst if required
    if gallery_conf["download_all_examples"]:
        download_fhindex = generate_zipfiles(
            gallery_dir_abs_path, app.builder.srcdir, gallery_conf
        )
        if sg_root_index:
            indexst += download_fhindex

    if sg_root_index:
        # Signature
        if app.config.sphinx_gallery_conf["show_signature"]:
            indexst += SPHX_GLR_SIG
        # Write index to file
        index_rst_new = os.path.join(gallery_dir_abs_path, "index.rst.new")
        with codecs.open(index_rst_new, "w", encoding="utf-8") as fhindex:
            fhindex.write(indexst)
        _replace_md5(index_rst_new, mode="t")


def _build_recommender(gallery_conf, gallery_dir_abs_path, subsecs):
    """Build recommender and write recommendations."""
    if gallery_conf["recommender"]["enable"]:
        try:
            import numpy as np  # noqa: F401
        except ImportError:
            raise ConfigError("gallery_conf['recommender'] requires numpy")

        recommender_params = copy.deepcopy(gallery_conf["recommender"])
        recommender_params.pop("enable")
        recommender_params.pop("rubric_header", None)
        recommender = ExampleRecommender(**recommender_params)

        gallery_py_files = []
        # root and subsection directories containing python examples
        gallery_directories = [gallery_dir_abs_path] + subsecs
        for current_dir in gallery_directories:
            src_dir = os.path.join(gallery_dir_abs_path, current_dir)
            # sort python files to have a deterministic input across call
            py_files = sorted(
                # NOTE we don't take account of `ignore_pattern` and ignore
                # ext in `example_extensions`
                [fname for fname in Path(src_dir).iterdir() if fname.suffix == ".py"],
                key=_get_callables(gallery_conf, "within_subsection_order", src_dir)[0],
            )
            gallery_py_files.append(
                [os.path.join(src_dir, fname) for fname in py_files]
            )
        # flatten the list of list
        gallery_py_files = list(chain.from_iterable(gallery_py_files))

        recommender.fit(gallery_py_files)
        for fname in gallery_py_files:
            _write_recommendations(recommender, fname, gallery_conf)


def _log_costs(costs, gallery_conf):
    """Log computation time."""
    logger.info("computation time summary:", color="white")
    lines, lens = _format_for_writing(
        costs, src_dir=gallery_conf["src_dir"], kind="console"
    )
    for name, t, m in lines:
        text = (f"    - {name}:   ").ljust(lens[0] + 10)
        if t is None:
            text += "(not run)"
            logger.info(text)
        else:
            t_float = float(t.split()[0])
            if t_float >= gallery_conf["min_reported_time"]:
                text += t.rjust(lens[1]) + "   " + m.rjust(lens[2])
                logger.info(text)


def generate_gallery_rst(app):
    """Generate the Main examples gallery reStructuredText.

    Fill Sphinx-Gallery configuration and scan example directories
    (up to one level depth of sub-directory) to generate example reST files.

    Iterate through each example directory and any of its sub-directories
    (creates sub-sections) that has a header/index file.
    Generate gallery example ReST files and `index.rst` file(s).

    If `nested_sections=True` we generate `index.rst` files for all
    sub-directories, which includes toctree linking to all sub-dir examples.
    The root example directory `index.rst` file will contain, in sequence,:

    * root gallery header then thumbnails,
    * toctree linking all examples in root gallery,
    * sub-section header followed by sub-section thumbnails, for all subsections,
    * a second final toctree, at the end of the file, linking to all sub-section
      index files.

    If `nested_sections=True` we generate a single `index.rst` file per
    example directory. It will contain headers for the root gallery and
    each sub-section, with each header followed by a toctree linking to
    every example in the root gallery/sub-section.
    """
    gallery_conf = app.config.sphinx_gallery_conf
    extra = ""
    if gallery_conf["parallel"]:
        extra = f" (with parallel={gallery_conf['parallel']})"
    logger.info(f"generating gallery{extra}...", color="white")

    seen_backrefs = set()

    costs = []
    workdirs = _prepare_sphx_glr_dirs(gallery_conf, app.builder.srcdir)

    # Check for duplicate filenames to make sure linking works as expected
    examples_dirs = [ex_dir for ex_dir, _ in workdirs]
    _collect_gallery_files(examples_dirs, gallery_conf, check_filenames=True)

    backrefs_all = {}

    for examples_dir, gallery_dir in workdirs:
        examples_dir_abs_path = os.path.join(app.builder.srcdir, examples_dir)
        gallery_dir_abs_path = os.path.join(app.builder.srcdir, gallery_dir)

        # Create example rst files for root gallery directory examples
        # (excl. sub-dir examples) and fetch gallery header for root index.rst
        (
            _,
            this_content,
            this_costs,
            this_toctree_items,
            backrefs_root,
        ) = generate_dir_rst(
            examples_dir_abs_path,
            gallery_dir_abs_path,
            gallery_conf,
            seen_backrefs,
            is_subsection=False,
        )

        _combine_backreferences(backrefs_all, backrefs_root)

        # `this_content` is None when user provides own index.rst
        sg_root_index = this_content is not None
        costs += this_costs
        write_computation_times(gallery_conf, gallery_dir_abs_path, this_costs)

        # `indexst` variable must exist, as passed to `_finish_index_rst`
        indexst = ""
        # Create root gallery index.rst
        if sg_root_index:
            # :orphan: to suppress "not included in TOCTREE" sphinx warnings
            indexst = ":orphan:\n\n" + this_content
            # Write toctree with gallery items from gallery root folder
            if len(this_toctree_items) > 0:
                this_toctree = _format_toctree(this_toctree_items)
                indexst += this_toctree

        # list all paths to subsection index files in this array
        subsection_index_files = []
        subsecs = get_subsections(
            app.builder.srcdir,
            examples_dir_abs_path,
            gallery_conf,
            check_for_header=sg_root_index,
        )
        for subsection in subsecs:
            src_dir = os.path.join(examples_dir_abs_path, subsection)
            target_dir = os.path.join(gallery_dir_abs_path, subsection)
            subsection_index_files.append(
                "/".join(["", gallery_dir, subsection, "index.rst"]).replace(
                    os.sep, "/"
                )  # fwd slashes needed in rst
            )

            (
                subsection_index_path,
                subsection_index_content,
                subsection_costs,
                subsection_toctree_filenames,
                backrefs_subsec,
            ) = generate_dir_rst(src_dir, target_dir, gallery_conf, seen_backrefs)

            _combine_backreferences(backrefs_all, backrefs_subsec)

            has_subsection_header = False
            if subsection_index_content:
                # Filter tags to prevent tag duplication across the documentation
                indexst += _filter_tags(subsection_index_content)
                has_subsection_header = True

            # Write subsection toctree containing all filenames, if req.
            if (
                sg_root_index
                and not gallery_conf["nested_sections"]
                and len(subsection_toctree_filenames) > 0
            ):
                subsection_index_toctree = _format_toctree(subsection_toctree_filenames)
                indexst += subsection_index_toctree
            # Otherwise, a new subsection index.rst.new file should
            # have been created and it needs to be parsed
            elif has_subsection_header:
                _replace_md5(subsection_index_path, mode="t")

            costs += subsection_costs
            write_computation_times(gallery_conf, target_dir, subsection_costs)

        # Per gallery - items below run once per gallery
        # Finish index.rst and write to file
        _finish_index_rst(
            app,
            gallery_conf,
            indexst,
            sg_root_index,
            subsection_index_files,
            gallery_dir_abs_path,
        )
        # Build recommendation system
        _build_recommender(gallery_conf, gallery_dir_abs_path, subsecs)

    # Per project - items below run once only (for all galleries)
    # Write a single global sg_execution_times
    write_computation_times(gallery_conf, None, costs)

    # Write backreferences_all to file
    if gallery_conf["backreferences_dir"]:
        _write_json(
            Path(
                gallery_conf["src_dir"],
                gallery_conf["backreferences_dir"],
                "backreferences_all",
            ),
            backrefs_all,
        )

    if gallery_conf["show_api_usage"] is not False:
        _init_api_usage(app.builder.srcdir)
    _finalize_backreferences(seen_backrefs, gallery_conf)

    if gallery_conf["plot_gallery"]:
        _log_costs(costs, gallery_conf)
        # Also create a junit.xml file, useful e.g. on CircleCI
        write_junit_xml(gallery_conf, app.builder.outdir, costs)


SPHX_GLR_ORPHAN = """
:orphan:

.. _{0}:

"""

SPHX_GLR_COMP_TIMES = (
    SPHX_GLR_ORPHAN
    + """
Computation times
=================
"""
)


def _sec_to_readable(t):
    """Convert a number of seconds to a more readable representation."""
    # This will only work for < 1 day execution time
    # And we reserve 2 digits for minutes because presumably
    # there aren't many > 99 minute scripts, but occasionally some
    # > 9 minute ones
    t = datetime(1, 1, 1) + timedelta(seconds=t)
    t = "{:02d}:{:02d}.{:03d}".format(
        t.hour * 60 + t.minute, t.second, int(round(t.microsecond / 1000.0))
    )
    return t


def _cost_key(cost):
    """Cost sorting function."""
    # sort by descending computation time, descending memory, alphabetical name
    return (-cost["t"], -cost["mem"], cost["src_file"])


def _format_for_writing(costs, *, src_dir, kind="rst"):
    """Provide formatted computation summary text.

    Parameters
    ----------
    costs: List[Dict]
        List of dicts of computation costs and paths, see gen_rst.py for details.
    src_dir : pathlib.Path
        The Sphinx source directory.
    kind: 'rst', 'rst-full' or 'console', default='rst'
        Format for printing to 'console' or for writing `sg_execution_times.rst' ('rst'
        for single galleries and 'rst-full' for all galleries).

    Returns
    -------
    lines: List[List[str]]
        Formatted computation text for each example, of format:
        [example_file, time_elapsed, memory_used]

    lens: List[int]
        Character length of each string in `lines`.
    """
    lines = list()
    for cost in sorted(costs, key=_cost_key):
        src_file = cost["src_file"]
        rel_path = os.path.relpath(src_file, src_dir)
        if kind in ("rst", "rst-full"):  # like in sg_execution_times
            target_dir_clean = os.path.relpath(cost["target_dir"], src_dir).replace(
                os.sep, "_"
            )
            paren = rel_path if kind == "rst-full" else os.path.basename(src_file)
            name = ":ref:`sphx_glr_{0}_{1}` (``{2}``)".format(
                target_dir_clean, os.path.basename(src_file), paren
            )
            t = _sec_to_readable(cost["t"])
        else:  # like in generate_gallery
            assert kind == "console"
            name = rel_path
            t = f"{cost['t']:0.2f} sec"
        m = f"{cost['mem']:.1f} MB"
        lines.append([name, t, m])
    lens = [max(x) for x in zip(*[[len(item) for item in cost] for cost in lines])]
    return lines, lens


def write_computation_times(gallery_conf, target_dir, costs):
    """Write computation times to `sg_execution_times.rst`.

    Parameters
    ----------
    gallery_conf : Dict[str, Any]
        Sphinx-Gallery configuration dictionary.
    target_dir : str | None
        Path to directory where example python source file are.
    costs: List[Dict]
        List of dicts of computation costs and paths, see gen_rst.py for details.
    """
    if not gallery_conf["write_computation_times"]:
        return
    total_time = sum(cost["t"] for cost in costs)
    if target_dir is None:  # all galleries together
        out_dir = gallery_conf["src_dir"]
        where = "all galleries"
        kind = "rst-full"
        ref_extra = ""
    else:  # a single gallery
        out_dir = target_dir
        where = os.path.relpath(target_dir, gallery_conf["src_dir"])
        kind = "rst"
        ref_extra = f"{where.replace(os.sep, '_')}_"
    new_ref = f"sphx_glr_{ref_extra}sg_execution_times"
    out_file = Path(out_dir) / "sg_execution_times.rst"
    if out_file.is_file() and total_time == 0:  # a re-run
        return
    with out_file.open("w", encoding="utf-8") as fid:
        fid.write(SPHX_GLR_COMP_TIMES.format(new_ref))
        fid.write(
            f"**{_sec_to_readable(total_time)}** total execution time for "
            f"{len(costs)} file{'s' if len(costs) != 1 else ''} **from {where}**:\n\n"
        )
        lines, lens = _format_for_writing(
            costs,
            src_dir=gallery_conf["src_dir"],
            kind=kind,
        )
        del costs
        # https://datatables.net/examples/styling/bootstrap5.html
        fid.write(  # put it in a container to make the scoped style work
            """\
.. container::

  .. raw:: html

    <style scoped>
    <link href="https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/5.3.0/css/bootstrap.min.css" rel="stylesheet" />
    <link href="https://cdn.datatables.net/1.13.6/css/dataTables.bootstrap5.min.css" rel="stylesheet" />
    </style>
    <script src="https://code.jquery.com/jquery-3.7.0.js"></script>
    <script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>
    <script src="https://cdn.datatables.net/1.13.6/js/dataTables.bootstrap5.min.js"></script>
    <script type="text/javascript" class="init">
    $(document).ready( function () {
        $('table.sg-datatable').DataTable({order: [[1, 'desc']]});
    } );
    </script>

  .. list-table::
   :header-rows: 1
   :class: table table-striped sg-datatable

   * - Example
     - Time
     - Mem (MB)
"""  # noqa: E501
        )
        # Need at least one entry or Sphinx complains
        for ex, t, mb in lines or [["N/A", "N/A", "N/A"]]:
            fid.write(
                f"""\
   * - {ex}
     - {t}
     - {mb.rsplit(maxsplit=1)[0]}
"""
            )  # remove the "MB" from the right


def write_api_entries(app, what, name, obj, options, lines):
    """Write api entries to `_sg_api_entries` configuration.

    To connect to `autodoc-process-docstring` event.

    Parameters
    ----------
    app :
        The Sphinx application object.
    what: str
        The type of the object which the docstring belongs to. One of
        "module", "class", "exception", "function", "method", "attribute".
    name :
        The fully qualified name of the object.
    obj :
        The object itself.
    options :
        The options given to the directive: an object with attributes inherited_members,
        undoc_members, show_inheritance and no-index that are true if the flag option of
        same name was given to the auto directive.
    lines :
        The lines of the docstring, see above.
    """
    if app.config.sphinx_gallery_conf["show_api_usage"] is False:
        return
    if "_sg_api_entries" not in app.config.sphinx_gallery_conf:
        app.config.sphinx_gallery_conf["_sg_api_entries"] = dict()
    if what not in app.config.sphinx_gallery_conf["_sg_api_entries"]:
        app.config.sphinx_gallery_conf["_sg_api_entries"][what] = set()
    app.config.sphinx_gallery_conf["_sg_api_entries"][what].add(name)


def _init_api_usage(gallery_dir):
    with codecs.open(
        os.path.join(gallery_dir, "sg_api_usage.rst"), "w", encoding="utf-8"
    ):
        pass


# Colors from https://personal.sron.nl/~pault/data/colourschemes.pdf
# 3 Diverging Colour Schemes, Figure 12, plus alpha=AA
API_COLORS = dict(
    edge="#00000080",  # gray (by alpha)
    okay="#98CAE180",  # blue
    bad_1="#FEDA8B80",  # yellow
    bad_2="#F67E4B80",  # orange
    bad_3="#A5002680",  # red
)


def _make_graph(fname, entries, gallery_conf):
    """Make a graph of unused and used API entries.

    The used API entries themselves are documented in the list, so
    for the graph, we'll focus on the number of unused API entries
    per modules. Modules with lots of unused entries (11+) will be colored
    red, those with less (6+) will be colored orange, those with only a few
    (1-5) will be colored yellow and those with no unused entries will be
    colored blue.

    The API entries that are used are shown with one graph per module.
    That way you can see the examples that each API entry is used in
    for that module (if this was done for the whole project at once,
    the graph would get too large very large quickly).

    Parameters
    ----------
    fname: str
        Path to '*sg_api_unused.dot' file.
    entries: Dict[str, List] or List[str]
        Used (List) or unused (Dict) API entries.
    gallery_conf : Dict[str, Any]
        Sphinx-Gallery configuration dictionary.
    """
    import graphviz

    dg = graphviz.Digraph(
        filename=fname,
        graph_attr={
            "overlap": "scale",
            "pad": "0.5",
        },
        node_attr={
            "color": API_COLORS["okay"],
            "style": "filled",
            "fontsize": "20",
            "shape": "box",
            "fontname": "Open Sans,Arial",
        },
    )

    if isinstance(entries, list):
        connections = set()
        lut = dict()  # look up table for connections so they don't repeat
        structs = [entry.split(".") for entry in entries]
        for struct in sorted(structs, key=len):
            for level in range(len(struct) - 2):
                if (struct[level], struct[level + 1]) in connections:
                    continue
                connections.add((struct[level], struct[level + 1]))
                node_from = (
                    lut[struct[level]] if struct[level] in lut else struct[level]
                )
                dg.node(node_from)
                node_to = struct[level + 1]
                node_kwargs = dict()
                # count, don't show leaves
                if len(struct) - 3 == level:
                    leaf_count = 0
                    for struct2 in structs:
                        # find structures of the same length as struct
                        if len(struct2) != level + 3:
                            continue
                        # find structures with two entries before
                        # the leaf that are the same as struct
                        if all(
                            [
                                struct2[level2] == struct[level2]
                                for level2 in range(level + 2)
                            ]
                        ):
                            leaf_count += 1
                    node_to += f" ({leaf_count})"
                    lut[struct[level + 1]] = node_to
                    if leaf_count > 10:
                        color_key = "bad_3"
                    elif leaf_count > 5:
                        color_key = "bad_2"
                    else:
                        color_key = "bad_1"
                    node_kwargs["color"] = API_COLORS[color_key]
                dg.node(node_to, **node_kwargs)
                dg.edge(node_from, node_to, color=API_COLORS["edge"])
        # add modules with all API entries
        for module in gallery_conf["_sg_api_entries"].get("module", []):
            struct = module.split(".")
            for i in range(len(struct) - 1):
                if struct[i + 1] not in lut:
                    dg.edge(struct[i], struct[i + 1])
    else:
        assert isinstance(entries, dict)
        for entry, refs in entries.items():
            dg.node(entry)
            for ref in refs:
                dg.node(ref, color=API_COLORS["bad_1"])
                dg.edge(entry, ref, color=API_COLORS["edge"])
    dg.save()


def write_api_entry_usage(app, docname, source):
    """Write an html page describing which API entries are used and unused.

    To document and graph only those API entries that are used by
    autodoc, we have to wait for autodoc to finish and hook into the
    ``source-read`` event. This intercepts the text from the rst such
    that it can be modified. Since, we only touched an empty file,
    we have to add 1) a list of all the API entries that are unused
    and a graph of the number of unused API entries per module and 2)
    a list of API entries that are used in examples, each with a sub-list
    of which examples that API entry is used in, and a graph that
    connects all of the API entries in a module to the examples
    that they are used in.

    Parameters
    ----------
    app :
        The Sphinx application object.
    docname :
        Docname of the document currently being parsed.
    source :
        List whose single element is the contents of the source file
    """
    docname = docname or ""  # can be None on Sphinx 7.2
    if docname != "sg_api_usage":
        return
    gallery_conf = app.config.sphinx_gallery_conf
    if gallery_conf["show_api_usage"] is False:
        return
    # since this is done at the gallery directory level (as opposed
    # to in a gallery directory, e.g. auto_examples), it runs last
    # which means that all the api entries will be in gallery_conf

    # Always write at least the title
    source[0] = SPHX_GLR_ORPHAN.format("sphx_glr_sg_api_usage")
    title = "Unused API Entries"
    source[0] += title + "\n" + "^" * len(title) + "\n\n"
    if (
        "_sg_api_entries" not in gallery_conf
        or gallery_conf["backreferences_dir"] is None
    ):
        source[0] += "No API entries found, not computed.\n\n"
        return
    backreferences_dir = os.path.join(
        gallery_conf["src_dir"], gallery_conf["backreferences_dir"]
    )

    example_files = set.union(
        *[
            gallery_conf["_sg_api_entries"][obj_type]
            for obj_type in ("class", "method", "function")
            if obj_type in gallery_conf["_sg_api_entries"]
        ]
    )

    if len(example_files) == 0:
        source[0] += "No examples run, not computed.\n\n"
        return

    def get_entry_type(entry):
        if entry in gallery_conf["_sg_api_entries"].get("class", []):
            return "class"
        elif entry in gallery_conf["_sg_api_entries"].get("method", []):
            return "meth"
        else:
            assert entry in gallery_conf["_sg_api_entries"]["function"]
            return "func"

    # find used and unused API entries
    unused_api_entries = list()
    used_api_entries = dict()
    backreferences_all = _read_json(Path(backreferences_dir, "backreferences_all.json"))
    src_dir = gallery_conf["src_dir"]
    for entry in example_files:
        # don't include built-in methods etc.
        if re.match(gallery_conf["api_usage_ignore"], entry) is not None:
            continue
        # check if backreferences empty
        backref_entry = backreferences_all.get(entry, None)
        if backref_entry is None:
            unused_api_entries.append(entry)
        else:
            used_api_entries[entry] = list()
            for br in backref_entry:
                # br[2] = abs path to target directory
                example_path = Path(br[2], br[0]).relative_to(src_dir)
                ref_name = str(example_path).replace(os.sep, "_")
                used_api_entries[entry].append(f"sphx_glr_{ref_name}")

    for entry in sorted(unused_api_entries):
        source[0] += f"- :{get_entry_type(entry)}:`{entry}`\n"
    source[0] += "\n\n"

    has_graphviz = _has_graphviz()
    if has_graphviz and unused_api_entries:
        source[0] += (
            ".. graphviz:: ./sg_api_unused.dot\n"
            "    :alt: API unused entries graph\n"
            "    :layout: neato\n\n"
        )

    used_count = len(used_api_entries)
    total_count = used_count + len(unused_api_entries)
    used_percentage = used_count / max(total_count, 1)  # avoid div by zero
    source[0] += (
        "\nAPI entries used: "
        f"{round(used_percentage * 100, 2)}% "
        f"({used_count}/{total_count})\n\n"
    )

    if has_graphviz and unused_api_entries:
        _make_graph(
            os.path.join(app.builder.srcdir, "sg_api_unused.dot"),
            unused_api_entries,
            gallery_conf,
        )

    if gallery_conf["show_api_usage"] is True and used_api_entries:
        title = "Used API Entries"
        source[0] += title + "\n" + "^" * len(title) + "\n\n"
        for entry in sorted(used_api_entries):
            source[0] += f"- :{get_entry_type(entry)}:`{entry}`\n\n"
            for ref in used_api_entries[entry]:
                source[0] += f"  - :ref:`{ref}`\n"
            source[0] += "\n\n"

        if has_graphviz:
            used_modules = {entry.split(".")[0] for entry in used_api_entries}
            for module in sorted(used_modules):
                source[0] += (
                    f"{module}\n" + "^" * len(module) + "\n\n"
                    f".. graphviz:: ./{module}_sg_api_used.dot\n"
                    f"    :alt: {module} usage graph\n"
                    "    :layout: neato\n\n"
                )

            for module in used_modules:
                logger.info("Making API usage graph for %s", module)
                # select and format entries for this module
                entries = dict()
                for entry, ref in used_api_entries.items():
                    if entry.split(".")[0] == module:
                        entry = entry.replace("sphx_glr_", "")
                        # remove prefix
                        for target_dir in gallery_conf["gallery_dirs"]:
                            if entry.startswith(target_dir):
                                entry = entry[len(target_dir) + 1 :]
                _make_graph(
                    os.path.join(app.builder.srcdir, f"{module}_sg_api_used.dot"),
                    entries,
                    gallery_conf,
                )


def clean_api_usage_files(app, exception):
    """Remove api usage .dot files.

    To connect to 'build-finished' event.
    """
    if os.path.isfile(os.path.join(app.builder.srcdir, "sg_api_usage.rst")):
        os.remove(os.path.join(app.builder.srcdir, "sg_api_usage.rst"))
    if os.path.isfile(os.path.join(app.builder.srcdir, "sg_api_unused.dot")):
        os.remove(os.path.join(app.builder.srcdir, "sg_api_unused.dot"))
    for file in os.listdir(app.builder.srcdir):
        if "sg_api_used.dot" in file:
            os.remove(os.path.join(app.builder.srcdir, file))


def write_junit_xml(gallery_conf, target_dir, costs):
    """Write JUnit XML file of example run times, successes, and failures.

    Parameters
    ----------
    gallery_conf : Dict[str, Any]
        Sphinx-Gallery configuration dictionary.
    target_dir : Union[str, pathlib.Path]
        Build directory.
    costs: List[Tuple[Tuple[float], str]]
        List of dicts of computation costs and paths, see gen_rst.py for details.
    """
    if not gallery_conf["junit"] or not gallery_conf["plot_gallery"]:
        return
    failing_as_expected, failing_unexpectedly, passing_unexpectedly = _parse_failures(
        gallery_conf
    )
    n_tests = 0
    n_failures = 0
    n_skips = 0
    elapsed = 0.0
    src_dir = gallery_conf["src_dir"]
    output = ""
    for cost in costs:
        t, fname = cost["t"], cost["src_file"]
        if not any(
            fname in x
            for x in (
                gallery_conf["passing_examples"],
                failing_unexpectedly,
                failing_as_expected,
                passing_unexpectedly,
            )
        ):
            continue  # not subselected by our regex
        title = gallery_conf["titles"][fname]
        output += (
            '<testcase classname={!s} file={!s} line="1" name={!s} time="{!r}">'.format(
                quoteattr(os.path.splitext(os.path.basename(fname))[0]),
                quoteattr(os.path.relpath(fname, src_dir)),
                quoteattr(title),
                t,
            )
        )
        if fname in failing_as_expected:
            output += '<skipped message="expected example failure"></skipped>'
            n_skips += 1
        elif fname in failing_unexpectedly or fname in passing_unexpectedly:
            if fname in failing_unexpectedly:
                traceback = gallery_conf["failing_examples"][fname]
            else:  # fname in passing_unexpectedly
                traceback = "Passed even though it was marked to fail"
            n_failures += 1
            output += "<failure message={!s}>{!s}</failure>".format(
                quoteattr(traceback.splitlines()[-1].strip()), escape(traceback)
            )
        output += "</testcase>"
        n_tests += 1
        elapsed += t
    output += "</testsuite>"
    output = (
        '<?xml version="1.0" encoding="utf-8"?>'
        '<testsuite errors="0" failures="{}" name="sphinx-gallery" '
        'skipped="{}" tests="{}" time="{}">'.format(
            n_failures, n_skips, n_tests, elapsed
        )
    ) + output
    # Actually write it
    fname = os.path.normpath(os.path.join(target_dir, gallery_conf["junit"]))
    junit_dir = os.path.dirname(fname)
    os.makedirs(junit_dir, exist_ok=True)
    with codecs.open(fname, "w", encoding="utf-8") as fid:
        fid.write(output)


def touch_empty_backreferences(app, what, name, obj, options, lines):
    """Generate empty back-reference example files.

    This avoids inclusion errors/warnings if there are no gallery
    examples for a class / module that is being parsed by autodoc.
    """
    if not bool(app.config.sphinx_gallery_conf["backreferences_dir"]):
        return

    examples_path = os.path.join(
        app.srcdir,
        app.config.sphinx_gallery_conf["backreferences_dir"],
        f"{name}.examples",
    )

    if not os.path.exists(examples_path):
        # touch file
        open(examples_path, "w").close()


def _expected_failing_examples(gallery_conf):
    return {
        os.path.normpath(os.path.join(gallery_conf["src_dir"], path))
        for path in gallery_conf["expected_failing_examples"]
    }


def _parse_failures(gallery_conf):
    """Split the failures."""
    failing_examples = set(gallery_conf["failing_examples"])
    expected_failing_examples = _expected_failing_examples(gallery_conf)
    failing_as_expected = failing_examples.intersection(expected_failing_examples)
    failing_unexpectedly = failing_examples.difference(expected_failing_examples)
    passing_unexpectedly = expected_failing_examples.difference(failing_examples)
    # filter from examples actually run
    passing_unexpectedly = set(
        src_file
        for src_file in passing_unexpectedly
        if re.search(gallery_conf["filename_pattern"], src_file)
    )
    return failing_as_expected, failing_unexpectedly, passing_unexpectedly


def summarize_failing_examples(app, exception):
    """Collects the list of falling examples and prints them with a traceback.

    Raises ValueError if there where failing examples.
    """
    if exception is not None:
        return

    # Under no-plot Examples are not run so nothing to summarize
    if not app.config.sphinx_gallery_conf["plot_gallery"]:
        logger.info(
            'Sphinx-Gallery gallery_conf["plot_gallery"] was '
            "False, so no examples were executed.",
            color="brown",
        )
        return

    gallery_conf = app.config.sphinx_gallery_conf
    failing_as_expected, failing_unexpectedly, passing_unexpectedly = _parse_failures(
        gallery_conf
    )

    idt = "    "
    if failing_as_expected:
        logger.info(
            bold(blue(f"Examples failing as expected ({len(failing_as_expected)}):"))
        )
        for fail_example in failing_as_expected:
            path = os.path.relpath(fail_example, gallery_conf["src_dir"])
            logger.info(
                f"{bold(blue(path))} failed leaving traceback:\n\n"
                f"{indent(gallery_conf['failing_examples'][fail_example], idt)}"
            )

    fail_msgs = []
    if failing_unexpectedly:
        fail_msgs.append(
            bold(red(f"Unexpected failing examples ({len(failing_unexpectedly)}):\n"))
        )
        for fail_example in failing_unexpectedly:
            path = os.path.relpath(fail_example, gallery_conf["src_dir"])
            fail_msgs.append(
                f"    {bold(red(path))} failed leaving traceback:\n\n"
                f"{indent(gallery_conf['failing_examples'][fail_example], idt)}"
            )

    if passing_unexpectedly:
        paths = [
            os.path.relpath(p, gallery_conf["src_dir"]) for p in passing_unexpectedly
        ]
        fail_msgs.append(
            bold(red(f"Examples expected to fail, but not failing ({len(paths)}):\n\n"))
            + red("\n".join(indent(p, idt) for p in paths))
            + "\n\nPlease remove these examples from "
            + "sphinx_gallery_conf['expected_failing_examples'] "
            + "in your conf.py file."
        )

    # standard message
    n_good = len(gallery_conf["passing_examples"])
    n_tot = len(gallery_conf["failing_examples"]) + n_good
    n_stale = len(gallery_conf["stale_examples"])
    logger.info(
        "\nSphinx-Gallery successfully executed %d out of %d "
        "file%s subselected by:\n\n"
        '    gallery_conf["filename_pattern"] = %r\n'
        '    gallery_conf["ignore_pattern"]   = %r\n'
        "\nafter excluding %d file%s that had previously been run "
        "(based on MD5).\n",
        n_good,
        n_tot,
        "s" if n_tot != 1 else "",
        gallery_conf["filename_pattern"],
        gallery_conf["ignore_pattern"],
        n_stale,
        "s" if n_stale != 1 else "",
        color="brown",
    )

    if fail_msgs:
        fail_message = bold(
            purple(
                "Here is a summary of the problems encountered "
                "when running the examples:\n\n"
                + "\n".join(fail_msgs)
                + "\n"
                + "-" * 79
            )
        )
        if gallery_conf["only_warn_on_example_error"]:
            logger.warning(fail_message)
        else:
            raise ExtensionError(fail_message)


def get_default_config_value(key):
    """Get default configuration function."""

    def default_getter(conf):
        return conf["sphinx_gallery_conf"].get(key, DEFAULT_GALLERY_CONF[key])

    return default_getter


def fill_gallery_conf_defaults(app, config, check_keys=True):
    """Check the sphinx-gallery config and set its defaults.

    This is called early at config-inited, so that all the rest of the code can
    do things like ``sphinx_gallery_conf['binder']['use_jupyter_lab']``, even
    if the keys have not been set explicitly in conf.py.
    """
    new_sphinx_gallery_conf = _fill_gallery_conf_defaults(
        config.sphinx_gallery_conf, app=app, check_keys=check_keys
    )
    config.sphinx_gallery_conf = new_sphinx_gallery_conf
    config.html_static_path.append(glr_path_static())

    # This must be done in "config-inited" for Sphinx 5 and 6 because otherwise the
    # `templates_path` will be overwritten and Sphinx will not be able to discover our
    # component templates
    config.templates_path.append(str(Path(__file__).parent / "components"))


def update_gallery_conf_builder_inited(app):
    """Update the the sphinx-gallery config at builder-inited."""
    plot_gallery = _bool_eval(app.builder.config.plot_gallery)
    src_dir = app.builder.srcdir
    abort_on_example_error = _bool_eval(app.builder.config.abort_on_example_error)
    _update_gallery_conf_builder_inited(
        app.config.sphinx_gallery_conf,
        src_dir,
        plot_gallery=plot_gallery,
        abort_on_example_error=abort_on_example_error,
        builder_name=app.builder.name,
    )


def setup_template_link_getters(app, pagename, templatename, context, doctree):
    """Set up the getters for download and launcher links.

    The getters are added to the sphinx context so as to be used in templates.
    """

    def _find_containers_with_class(class_name):
        if doctree is None:
            return iter([])

        return doctree.findall(
            lambda x: isinstance(x, nodes.container)
            and class_name in x.attributes.get("classes", [])
        )

    def get_download_links():
        """Get the download links for the example.

        This function relies on `_find_containers_with_class` which in turn relies on
        `doctree` provided in the `html-page-context` event. This function will then
        be added to the context of the page and can be accessed in the template.

        This returns a dictionary with keys in ["python", "jupyter", "zip"], depending
        on their availability. The values contain:
        - link: The relative path to the download file
        - label: The "Download {label}" text
        - title: The title to show when hovering over the link
        """
        links = {}
        for key, label in [
            ("python", "source code"),
            ("jupyter", "Jupyter notebook"),
            ("zip", "zipped"),
        ]:
            containers = _find_containers_with_class(f"sphx-glr-download-{key}")
            if container := next(containers, None):
                attrs = container.children[0].children[0].attributes
                if link := attrs.get("filename"):
                    links[key] = {
                        "link": f"_downloads/{link}",
                        "label": label,
                        "title": attrs.get("reftarget"),
                    }
        return links

    def get_launcher_links():
        """Get the launcher links for the example.

        This function relies on `_find_containers_with_class` which in turn relies on
        `doctree` provided in the `html-page-context` event. This function will then
        be added to the context of the page and can be accessed in the template.

        This returns a dictionary with keys in ["lite", "binder"], depending on their
        availability. The values contain:
        - link: The URL to Binder or JupyterLite link of the example
        - img_alt: The alt text for the link badge
        - img_src: The source URL of the link badge
        """
        links = {}
        for key in ["lite", "binder"]:
            containers = _find_containers_with_class(f"{key}-badge")
            if container := next(containers, None):
                anchor = container.children[0]
                image = anchor.children[0]
                if link := anchor.attributes.get("refuri"):
                    links[key] = {
                        "link": link,
                        "img_alt": image.attributes.get("alt"),
                        "img_src": image.attributes.get("uri"),
                    }
        return links

    context["get_download_links"] = get_download_links
    context["get_launcher_links"] = get_launcher_links


def setup(app):
    """Setup Sphinx-Gallery sphinx extension."""
    app.add_config_value("sphinx_gallery_conf", DEFAULT_GALLERY_CONF, "html")
    for key in ["plot_gallery", "abort_on_example_error"]:
        app.add_config_value(key, get_default_config_value(key), "html")

    # Early filling of sphinx_gallery_conf defaults at config-inited
    app.connect("config-inited", fill_gallery_conf_defaults, priority=10)
    # set small priority value, so that pre_configure_jupyterlite_sphinx is
    # called before jupyterlite_sphinx config-inited
    app.connect("config-inited", pre_configure_jupyterlite_sphinx, priority=100)
    # set high priority value, so that post_configure_jupyterlite_sphinx is
    # called after jupyterlite_sphinx config-inited
    app.connect("config-inited", post_configure_jupyterlite_sphinx, priority=900)

    if "sphinx.ext.autodoc" in app.extensions:
        app.connect("autodoc-process-docstring", touch_empty_backreferences)
        app.connect("autodoc-process-docstring", write_api_entries)
        app.connect("source-read", write_api_entry_usage)

    # Add the custom directive
    app.add_directive("minigallery", MiniGallery)
    app.add_directive("image-sg", ImageSg)

    imagesg_addnode(app)

    # Early update of sphinx_gallery_conf at builder-inited
    app.connect("builder-inited", update_gallery_conf_builder_inited, priority=10)
    app.connect("builder-inited", generate_gallery_rst)
    app.connect("build-finished", copy_binder_files)
    app.connect("build-finished", create_jupyterlite_contents)

    app.connect("build-finished", summarize_failing_examples)
    app.connect("build-finished", embed_code_links)
    app.connect("build-finished", clean_api_usage_files)

    app.connect("html-page-context", setup_template_link_getters)

    metadata = {
        "parallel_read_safe": True,
        "parallel_write_safe": True,
        "version": _sg_version,
    }
    return metadata


def setup_module():
    """Hack to stop nosetests running setup() above."""
    pass
