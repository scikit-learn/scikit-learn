# Author: Óscar Nájera
# License: 3-clause BSD
r"""Test Sphinx-Gallery gallery generation."""

import inspect
import json
import os
import re
from pathlib import Path

import pytest
from sphinx.application import Sphinx
from sphinx.config import is_serializable
from sphinx.errors import ConfigError, ExtensionError, SphinxWarning

from sphinx_gallery.gen_gallery import (
    _fill_gallery_conf_defaults,
    fill_gallery_conf_defaults,
    write_api_entry_usage,
    write_computation_times,
)
from sphinx_gallery.interactive_example import create_jupyterlite_contents
from sphinx_gallery.utils import (
    _collect_gallery_files,
    _escape_ansi,
    check_duplicate_filenames,
    check_spaces_in_filenames,
)

MINIMAL_HEADER = """
'''
Title
-----
Description.
'''

"""


def test_bad_config():
    """Test that bad config values are caught."""
    sphinx_gallery_conf = dict(example_dir="")
    with pytest.raises(
        ConfigError, match="example_dir.*did you mean 'examples_dirs'?.*"
    ):
        _fill_gallery_conf_defaults(sphinx_gallery_conf)
    sphinx_gallery_conf = dict(n_subsection_order="")
    with pytest.raises(
        ConfigError, match=r"did you mean one of \['subsection_order', 'within_.*"
    ):
        _fill_gallery_conf_defaults(sphinx_gallery_conf)
    sphinx_gallery_conf = dict(within_subsection_order="sphinx_gallery.a.b.Key")
    with pytest.raises(ConfigError, match=r"Unknown string option.*when importing"):
        _fill_gallery_conf_defaults(sphinx_gallery_conf)
    sphinx_gallery_conf = dict(within_subsection_order=1.0)
    with pytest.raises(ConfigError, match="must be callable"):
        _fill_gallery_conf_defaults(sphinx_gallery_conf)
    sphinx_gallery_conf = dict(minigallery_sort_order=int)
    with pytest.raises(ConfigError, match="Got class rather than callable instance"):
        _fill_gallery_conf_defaults(sphinx_gallery_conf)


def test_default_config(sphinx_app_wrapper):
    """Test default Sphinx-Gallery config loaded when extension added to Sphinx."""
    sphinx_app = sphinx_app_wrapper.create_sphinx_app()
    cfg = sphinx_app.config
    assert cfg.project == "Sphinx-Gallery <Tests>"
    # no duplicate values allowed The config is present already
    with pytest.raises(ExtensionError) as excinfo:
        sphinx_app.add_config_value("sphinx_gallery_conf", "x", True)
    assert "already present" in str(excinfo.value)


def test_serializable(sphinx_app_wrapper):
    """Test that the default config is serializable."""
    bad = list()
    for key, val in _fill_gallery_conf_defaults({}).items():
        if not is_serializable(val):
            bad.append(f"{repr(key)}: {repr(val)}")

    bad = "\n".join(bad)
    assert not bad, f"Non-serializable values found:\n{bad}"


@pytest.mark.add_conf(
    content="""
sphinx_gallery_conf = {
    'examples_dirs': 'src',
    'gallery_dirs': 'ex',
}"""
)
def test_no_warning_simple_config(sphinx_app_wrapper):
    """Testing that no warning is issued with a simple config.

    The simple config only specifies input (examples_dirs) and output
    (gallery_dirs) directories.
    """
    sphinx_app = sphinx_app_wrapper.create_sphinx_app()
    cfg = sphinx_app.config
    assert cfg.project == "Sphinx-Gallery <Tests>"
    build_warn = sphinx_app._warning.getvalue()
    assert build_warn == ""


@pytest.mark.add_conf(
    content="""
sphinx_gallery_conf = {
    'examples_dirs': ['src', 'excess'],
    'gallery_dirs': 'ex',
}"""
)
def test_unequal_examples_gallery_dirs(sphinx_app_wrapper):
    """Check warning when 'examples_dirs' and 'gallery_dirs' unequal lengths."""
    sphinx_app = sphinx_app_wrapper.create_sphinx_app()
    build_warn = sphinx_app._warning.getvalue()
    assert "'examples_dirs' and 'gallery_dirs' are of different lengths" in build_warn


# This changed from passing the ValueError directly to
# raising "sphinx.errors.ConfigError" with "threw an exception"
@pytest.mark.parametrize(
    "err_class, err_match",
    [
        pytest.param(
            ConfigError,
            "Unknown string option for reset_modules",
            id="Resetter unknown",
            marks=pytest.mark.add_conf(
                content="sphinx_gallery_conf={'reset_modules': ('f',)}"
            ),
        ),
        pytest.param(
            ConfigError,
            "reset_modules.* must be callable",
            id="Resetter not callable",
            marks=pytest.mark.add_conf(
                content="sphinx_gallery_conf={'reset_modules': (1.,),}"
            ),
        ),
    ],
)
def test_bad_reset(sphinx_app_wrapper, err_class, err_match):
    with pytest.raises(err_class, match=err_match):
        sphinx_app_wrapper.create_sphinx_app()


@pytest.mark.parametrize(
    "err_class, err_match",
    [
        pytest.param(
            ConfigError,
            "'reset_modules_order' config allowed types",
            id="Resetter unknown",
            marks=pytest.mark.add_conf(
                content=("sphinx_gallery_conf={'reset_modules_order': 1,}")
            ),
        ),
        pytest.param(
            ConfigError,
            "reset_modules_order must be in",
            id="reset_modules_order not valid",
            marks=pytest.mark.add_conf(
                content=("sphinx_gallery_conf={'reset_modules_order': 'invalid',}")
            ),
        ),
    ],
)
def test_bad_reset_modules_order(sphinx_app_wrapper, err_class, err_match):
    with pytest.raises(err_class, match=err_match):
        sphinx_app_wrapper.create_sphinx_app()


@pytest.mark.parametrize(
    "err_class, err_match",
    [
        pytest.param(
            ConfigError,
            "Unknown css",
            id="CSS str error",
            marks=pytest.mark.add_conf(content="sphinx_gallery_conf={'css': ('foo',)}"),
        ),
        pytest.param(
            ConfigError,
            "config allowed types:",
            id="CSS type error",
            marks=pytest.mark.add_conf(content="sphinx_gallery_conf={'css': 1.}"),
        ),
    ],
)
def test_bad_css(sphinx_app_wrapper, err_class, err_match):
    """Test 'css' configuration validation is correct."""
    with pytest.raises(err_class, match=err_match):
        sphinx_app_wrapper.create_sphinx_app()


def test_bad_api():
    """Test that we raise an error for bad API usage arguments."""
    sphinx_gallery_conf = dict(api_usage_ignore=("foo",))
    with pytest.raises(ConfigError, match="'api_usage_ignore' config allowed"):
        _fill_gallery_conf_defaults(sphinx_gallery_conf)
    sphinx_gallery_conf = dict(show_api_usage="foo")
    with pytest.raises(ConfigError, match='.*must be True, False or "unused".*'):
        _fill_gallery_conf_defaults(sphinx_gallery_conf)


@pytest.mark.add_conf(
    content="""
sphinx_gallery_conf = {
    'backreferences_dir': os.path.join('gen_modules', 'backreferences'),
    'examples_dirs': 'src',
    'gallery_dirs': 'ex',
}"""
)
def test_config_backreferences(sphinx_app_wrapper):
    """Test no warning is issued under the new configuration."""
    sphinx_app = sphinx_app_wrapper.create_sphinx_app()
    cfg = sphinx_app.config
    assert cfg.project == "Sphinx-Gallery <Tests>"
    assert cfg.sphinx_gallery_conf["backreferences_dir"] == os.path.join(
        "gen_modules", "backreferences"
    )
    build_warn = sphinx_app._warning.getvalue()
    assert build_warn == ""


def test_duplicate_files_warn(sphinx_app_wrapper):
    """Test for a warning when two files with the same filename exist."""
    sphinx_app = sphinx_app_wrapper.create_sphinx_app()

    files = ["./a/file1.py", "./a/file2.py", "a/file3.py", "./b/file1.py"]
    msg = (
        "Duplicate example file name(s) found. Having duplicate file names "
        "will break some links. List of files: {}"
    )
    m = "['./b/file1.py']"

    # No warning because no overlapping names
    check_duplicate_filenames(files[:-1])
    build_warn = sphinx_app._warning.getvalue()
    assert build_warn == ""

    # Warning because last file is named the same
    check_duplicate_filenames(files)
    build_warn = sphinx_app._warning.getvalue()
    assert msg.format(m) in build_warn


def test_spaces_in_files_warn(sphinx_app_wrapper):
    """Test for a exception when an example filename has a space in it."""
    sphinx_app = sphinx_app_wrapper.create_sphinx_app()

    files = ["./a/file1.py", "./a/file2.py", "./a/file 3.py"]
    msg = (
        "Example file name(s) with spaces found. Having spaces in "
        "file names will break some links. "
        "List of files: {}"
    )
    m = "['./a/file 3.py']"

    # No warning because no filename with space
    check_spaces_in_filenames(files[:-1])
    build_warn = sphinx_app._warning.getvalue()
    assert build_warn == ""

    # Warning because last file has space
    check_spaces_in_filenames(files)
    build_warn = sphinx_app._warning.getvalue()
    assert msg.format(m) in build_warn


def _check_order(sphinx_app, key, expected_order=None):
    """Iterates through sphx-glr-thumbcontainer divs and reads key from the tooltip.

    Used to test that these keys (in index.rst) appear in a specific order. The `.py`
    files include a source-of-truth line indicating the position in the sort order that
    the thumbnail for that file should occur, for the various built-in sorters.

    The regex below extracts info from that source-of-truth line to check that the
    thumbnail ordering was correct --- unless `expected_order` is provided, in which
    case an alternative regex extracts the digit-part of the filename from the crossref,
    and compares the order of crossrefs to `expected_order`.
    """
    index_fname = Path(sphinx_app.outdir, "..", "ex", "index.rst")
    order = list()
    if key is None:
        regex = r".*:ref:`sphx_glr_ex_(?:(?:first|second)-subsection_)?plot_(\d)\.py`"
        locator = ":ref:"
    else:
        regex = rf".*:{key}=(\d):.*"
        locator = "sphx-glr-thumbcontainer"
    with open(index_fname, "r", encoding="utf-8") as fid:
        for line in fid:
            if locator in line:
                order.append(re.match(regex, line).group(1))
    expected_order = expected_order or list("123456789")
    assert order == expected_order


# ======================================================================================
# SUBSECTION AND WITHIN-SUBSECTION SORTING TESTS
# ======================================================================================
# pre-construct the pytest params for the builtin sorters, as they're all very similar
_template_conf_for_builtin_sorters = """
sphinx_gallery_conf = {{
    "examples_dirs": "src",
    "gallery_dirs": "ex",
    "within_subsection_order": "{}",
}}"""
_builtin_sorters = dict(
    lines="NumberOfCodeLinesSortKey",  # the default
    filesize="FileSizeSortKey",
    filename="FileNameSortKey",
    title="ExampleTitleSortKey",
)
_params_for_testing_builtin_sorters = (
    pytest.param(
        kind,  # `sort_key` param
        None,  # `expected_order` (order is extracted from comments in test site pages)
        id=f"within_subsection_sort_by_{kind}",
        marks=pytest.mark.add_conf(
            content=_template_conf_for_builtin_sorters.format(sorter)
        ),
    )
    for kind, sorter in _builtin_sorters.items()
)
# now a more flexible config template, for changing subsection and/or within-subsection
# order (possibly at the same time):
_template_conf = """
{imports}
sphinx_gallery_conf = {{
    "examples_dirs": "src",
    "gallery_dirs": "ex",
    "subsection_order": {subsection_order},
    "within_subsection_order": {within_subsection_order},
}}"""
# fill in the template config with different test cases and wrap in an `add_conf` mark
_subsection_explicit_order = pytest.mark.add_conf(
    content=_template_conf.format(
        imports="from sphinx_gallery.sorting import ExplicitOrder",
        subsection_order='ExplicitOrder(["src/second-subsection", "src/first-subsection"])',
        within_subsection_order='"NumberOfCodeLinesSortKey"',  # this is the default
    )
)
_subsection_explicit_order_list = pytest.mark.add_conf(
    content=_template_conf.format(
        imports="",
        subsection_order='["src/second-subsection", "src/first-subsection"]',
        within_subsection_order="sphinx_gallery.sorting.NumberOfCodeLinesSortKey",
    )
)
_both_custom_fqn = pytest.mark.add_conf(
    content=_template_conf.format(
        imports="",
        subsection_order='"sphinx_gallery.utils._custom_subsection_sorter"',
        within_subsection_order='"sphinx_gallery.utils._custom_example_sorter"',
    )
)
_within_subsection_custom_fqn = pytest.mark.add_conf(
    content=_template_conf.format(
        imports="",
        subsection_order="None",  # the default, AKA, sort by foldername
        within_subsection_order='"sphinx_gallery.utils._custom_example_sorter"',
    )
)
_custom_func = """
from sphinx_gallery.sorting import FunctionSortKey

def custom_sorter(filename):
    ORDER = [
        "plot_3.py",
        "plot_2.py",
        "plot_1.py",
        "plot_6.py",
        "plot_5.py",
        "plot_4.py",
        "plot_9.py",
        "plot_7.py",
        "plot_8.py",
    ]
    return ORDER.index(filename)
"""
_within_subsection_custom_func = pytest.mark.add_conf(
    content=_template_conf.format(
        imports=_custom_func,
        subsection_order="None",  # the default, AKA, sort by foldername
        within_subsection_order="FunctionSortKey(custom_sorter)",
    )
)
_subsection_fqn_within_subsection_custom_func = pytest.mark.add_conf(
    content=_template_conf.format(
        imports=_custom_func,
        subsection_order='"sphinx_gallery.utils._custom_subsection_sorter"',
        within_subsection_order="FunctionSortKey(custom_sorter)",
    )
)


@pytest.mark.parametrize(
    "sort_key,expected_order",
    (
        pytest.param(
            "lines",
            # ↓ natural order except subsection 2 (789) precedes subsection 1 (456)
            list("123789456"),
            id="subsection_sort_by_ExplicitOrder",
            marks=_subsection_explicit_order,
        ),
        pytest.param(
            "lines",
            # ↓ natural order except subsection 2 (789) precedes subsection 1 (456)
            list("123789456"),
            id="subsection_sort_by_list",
            marks=_subsection_explicit_order_list,
        ),
        pytest.param(
            None,
            # ↓ order determined by `utils._CUSTOM_EXAMPLE_ORDER`
            list("132564879"),
            id="within_subsection_sort_by_custom_FQN",
            marks=_within_subsection_custom_fqn,
        ),
        pytest.param(
            None,
            # ↓ order set by `utils._CUSTOM_EXAMPLE_ORDER`, with middle (564) and
            # ↓ last (879) thirds swapped due to subsection reordering
            list("132879564"),
            id="subsection_and_within_subsection_both_sort_by_custom_FQN",
            marks=_both_custom_fqn,
        ),
        pytest.param(
            None,
            # ↓ order set by local `_custom_func` variable within this test file
            list("321654978"),
            id="within_subsection_sort_by_custom_func",
            marks=_within_subsection_custom_func,
        ),
        pytest.param(
            None,
            # ↓ order set by `_custom_func` in this test file, with middle (654) and
            # ↓ last (978) thirds swapped due to subsection reordering
            list("321978654"),
            id="subsection_sort_by_FQN_and_within_subsection_sort_by_custom_func",
            marks=_subsection_fqn_within_subsection_custom_func,
        ),
        *_params_for_testing_builtin_sorters,
    ),
)
def test_example_sorting(sphinx_app_wrapper, sort_key, expected_order):
    """Test sorting of examples with fully qualified name of a custom sort function."""
    sphinx_app = sphinx_app_wrapper.create_sphinx_app()
    _check_order(sphinx_app, sort_key, expected_order=expected_order)


def test_collect_gallery_files(tmpdir, gallery_conf):
    """Test that example files are collected properly."""
    rel_filepaths = [
        "examples/file1.py",
        "examples/test.rst",
        "examples/GALLERY_HEADER.rst",
        "examples/folder1/file1.py",
        "examples/folder1/file2.py",
        "examples/folder2/file1.py",
        "tutorials/folder1/file1.py",
        "tutorials/folder2/file1.py",
    ]

    abs_paths = [tmpdir.join(rp) for rp in rel_filepaths]
    for ap in abs_paths:
        ap.ensure()

    examples_path = tmpdir.join("examples")
    dirs = [examples_path.strpath]
    collected_files = set(
        _collect_gallery_files(dirs, gallery_conf, check_filenames=True)
    )
    expected_files = {
        ap.strpath for ap in abs_paths if re.search(r"examples.*\.py$", ap.strpath)
    }

    assert collected_files == expected_files

    tutorials_path = tmpdir.join("tutorials")
    dirs = [examples_path.strpath, tutorials_path.strpath]
    collected_files = set(
        _collect_gallery_files(dirs, gallery_conf, check_filenames=True)
    )
    expected_files = {
        ap.strpath for ap in abs_paths if re.search(r".*\.py$", ap.strpath)
    }

    assert collected_files == expected_files


def test_collect_gallery_files_ignore_pattern(tmpdir, gallery_conf):
    """Test that ignore pattern example files are not collected."""
    rel_filepaths = [
        "examples/file1.py",
        "examples/folder1/fileone.py",
        "examples/folder1/file2.py",
        "examples/folder2/fileone.py",
    ]

    abs_paths = [tmpdir.join(rp) for rp in rel_filepaths]
    for ap in abs_paths:
        ap.ensure()

    gallery_conf["ignore_pattern"] = r"one"
    examples_path = tmpdir.join("examples")
    dirs = [examples_path.strpath]
    collected_files = set(
        _collect_gallery_files(dirs, gallery_conf, check_filenames=True)
    )
    expected_files = {
        ap.strpath
        for ap in abs_paths
        if re.search(r"one", Path(ap.strpath).name) is None
    }

    assert collected_files == expected_files


@pytest.mark.add_conf(
    content=r"""
sphinx_gallery_conf = {
    'examples_dirs': 'src',
    'gallery_dirs': 'ex',
    'copyfile_regex': r'.*\.rst',
}"""
)
@pytest.mark.add_rst(file="own index.rst")
def test_own_index_first(sphinx_app_wrapper):
    """Test `generate_gallery_rst` works when own index gallery is first (and only)."""
    # Issue #1382
    sphinx_app_wrapper.build_sphinx_app()


@pytest.mark.add_conf(
    content="""
sphinx_gallery_conf = {
    'backreferences_dir' : os.path.join('modules', 'gen'),
    'examples_dirs': 'src',
    'gallery_dirs': ['ex'],
    'binder': {'binderhub_url': 'http://test1.com', 'org': 'org',
               'repo': 'repo', 'branch': 'branch',
               'notebooks_dir': 'ntbk_folder',
               'dependencies': 'requirements.txt'}
}"""
)
def test_binder_copy_files(sphinx_app_wrapper):
    """Test that notebooks are copied properly."""
    from sphinx_gallery.interactive_example import copy_binder_files

    sphinx_app = sphinx_app_wrapper.create_sphinx_app()
    gallery_conf = sphinx_app.config.sphinx_gallery_conf
    # Create requirements file
    with open(Path(sphinx_app.srcdir, "requirements.txt"), "w"):
        pass
    copy_binder_files(sphinx_app, None)

    for i_file in ["plot_1", "plot_2", "plot_3"]:
        assert Path(
            sphinx_app.outdir,
            "ntbk_folder",
            gallery_conf["gallery_dirs"][0],
            i_file + ".ipynb",
        ).exists()


@pytest.mark.add_conf(
    content="""
sphinx_gallery_conf = {
    'examples_dirs': 'src',
    'gallery_dirs': 'ex',
}"""
)
def test_failing_examples_raise_exception(sphinx_app_wrapper):
    example_dir = Path(sphinx_app_wrapper.srcdir, "src")
    bad_line = "print(f'{a[}')"  # never closed bracket inside print -> SyntaxError
    bad_code = f"""\
'''
Failing example
---------------
Should emit a syntax error in the second code block.
'''
1 + 2

# %%
# More

{bad_line}
"""
    bad_line_no = bad_code.split("\n").index(bad_line) + 1
    with open(Path(example_dir, "plot_3.py"), "w", encoding="utf-8") as fid:
        fid.write(bad_code)
    with pytest.raises(ExtensionError) as excinfo:
        sphinx_app_wrapper.build_sphinx_app()
    tb = str(excinfo.value)
    assert "Unexpected failing examples" in tb
    # Check traceback points to correct line (see #1301)
    assert f"line {bad_line_no}" in tb
    assert bad_line in tb


@pytest.mark.add_conf(
    content="""
sphinx_gallery_conf = {
    'examples_dirs': 'src',
    'gallery_dirs': 'ex',
    'filename_pattern': 'plot_1.py',
}"""
)
def test_expected_failing_examples_were_executed(sphinx_app_wrapper):
    """Testing that no exception is issued when broken example is not built.

    See #335 for more details.
    """
    sphinx_app_wrapper.build_sphinx_app()


@pytest.mark.add_conf(
    content="""
sphinx_gallery_conf = {
    'examples_dirs': 'src',
    'gallery_dirs': 'ex',
    'only_warn_on_example_error': True,
}"""
)
def test_only_warn_on_example_error(sphinx_app_wrapper):
    """Test behaviour of only_warn_on_example_error flag."""
    example_dir = Path(sphinx_app_wrapper.srcdir) / "src"
    with open(example_dir / "plot_3.py", "w", encoding="utf-8") as fid:
        fid.write(f"{MINIMAL_HEADER}raise ValueError")
    sphinx_app = sphinx_app_wrapper.build_sphinx_app()

    build_warn = _escape_ansi(sphinx_app._warning.getvalue())
    assert "plot_3.py unexpectedly failed to execute correctly" in build_warn
    assert "WARNING: Here is a summary of the problems" in build_warn


@pytest.mark.add_conf(
    content="""
sphinx_gallery_conf = {
    'examples_dirs': 'src',
    'gallery_dirs': 'ex',
    'only_warn_on_example_error': True,
}"""
)
def test_only_warn_on_example_error_sphinx_warning(sphinx_app_wrapper):
    """Test behaviour of only_warn_on_example_error flag."""
    # https://github.com/sphinx-doc/sphinx/pull/12743/files
    for key in ("warningiserror", "exception_on_warning"):
        if key in inspect.getfullargspec(Sphinx).args:
            sphinx_app_wrapper.kwargs[key] = True
    example_dir = Path(sphinx_app_wrapper.srcdir) / "src"
    with open(example_dir / "plot_3.py", "w", encoding="utf-8") as fid:
        fid.write(f"{MINIMAL_HEADER}raise ValueError")
    with pytest.raises(SphinxWarning) as excinfo:
        sphinx_app_wrapper.build_sphinx_app()
    exc = _escape_ansi(str(excinfo.value))
    assert "plot_3.py unexpectedly failed to execute" in exc


@pytest.mark.add_conf(
    content="""
sphinx_gallery_conf = {
    'examples_dirs': 'src',
    'gallery_dirs': 'ex',
    'expected_failing_examples' :['src/plot_2.py'],
}"""
)
def test_examples_not_expected_to_pass(sphinx_app_wrapper):
    with pytest.raises(ExtensionError) as excinfo:
        sphinx_app_wrapper.build_sphinx_app()
    exc = _escape_ansi(str(excinfo.value))
    assert "expected to fail, but not failing" in exc


@pytest.mark.add_conf(
    content="""
from sphinx_gallery.gen_rst import _sg_call_memory_noop

sphinx_gallery_conf = {
    'show_memory': _sg_call_memory_noop,
    'examples_dirs': 'src',
    'gallery_dirs': 'ex',
}"""
)
def test_show_memory_callable(sphinx_app_wrapper):
    sphinx_app = sphinx_app_wrapper.build_sphinx_app()
    status = sphinx_app._status.getvalue()
    assert "0.0 MB" in status, status


@pytest.mark.parametrize(
    "",
    [
        pytest.param(
            id="first notebook cell",
            marks=pytest.mark.add_conf(
                content="""sphinx_gallery_conf = {'first_notebook_cell': 2,}"""
            ),
        ),
        pytest.param(
            id="last notebook cell",
            marks=pytest.mark.add_conf(
                content="""sphinx_gallery_conf = {'last_notebook_cell': 2,}"""
            ),
        ),
    ],
)
def test_notebook_cell_config(sphinx_app_wrapper):
    """Tests that first and last cell configuration validated."""
    with pytest.raises(ConfigError):
        app = sphinx_app_wrapper.create_sphinx_app()
        fill_gallery_conf_defaults(app, app.config, check_keys=False)


@pytest.mark.add_conf(
    content="""
sphinx_gallery_conf = {
    'backreferences_dir': False,
}"""
)
def test_backreferences_dir_config(sphinx_app_wrapper):
    """Tests 'backreferences_dir' type checking."""
    with pytest.raises(ConfigError, match="'backreferences_dir' config allowed types"):
        app = sphinx_app_wrapper.create_sphinx_app()
        fill_gallery_conf_defaults(app, app.config, check_keys=False)


@pytest.mark.add_conf(
    content="""
sphinx_gallery_conf = {
    'examples_dirs': 'src',
    'gallery_dirs': 'ex',
}"""
)
def test_minigallery_no_backreferences_dir(sphinx_app_wrapper):
    """Check warning when no backreferences_dir set but minigallery directive used."""
    sphinx_app = sphinx_app_wrapper.build_sphinx_app()
    build_warn = sphinx_app._warning.getvalue()
    assert "'backreferences_dir' config is None, minigallery" in build_warn


@pytest.mark.add_conf(
    content="""
import pathlib

sphinx_gallery_conf = {
    'backreferences_dir': pathlib.Path('.'),
}"""
)
def test_backreferences_dir_pathlib_config(sphinx_app_wrapper):
    """Tests pathlib.Path does not raise exception."""
    app = sphinx_app_wrapper.create_sphinx_app()
    fill_gallery_conf_defaults(app, app.config, check_keys=False)


@pytest.mark.add_conf(
    content="""
sphinx_gallery_conf = {
    'examples_dirs': 'src',
    'gallery_dirs': 'ex',
}"""
)
@pytest.mark.add_rst(
    file="""
Header
======

.. minigallery:: index.rst
"""
)
def test_minigallery_not_in_examples_dirs(sphinx_app_wrapper):
    """Check error when minigallery directive's path input not in `examples_dirs`."""
    msg = "minigallery directive error: path input 'index.rst'"
    with pytest.raises(ExtensionError, match=msg):
        sphinx_app_wrapper.build_sphinx_app()


@pytest.mark.add_conf(
    content="""
sphinx_gallery_conf = {
    'examples_dirs': ['src', 'src/sub_folder/sub_sub_folder'],
    'gallery_dirs': ['ex', 'ex/sub_folder/sub_sub_folder'],
}"""
)
@pytest.mark.add_rst(
    file="""
Header
======

.. minigallery:: src/sub_folder/sub_sub_folder/plot_nested.py
"""
)
def test_minigallery_multi_match(sphinx_app_wrapper):
    """Check minigallery directive's path input resolution in nested `examples_dirs`.

    When a examples gallery is nested inside another examples gallery, path inputs
    from the nested gallery should resolve to the nested gallery.
    """
    sphinx_app = sphinx_app_wrapper.build_sphinx_app()
    minigallery_html = Path(sphinx_app.outdir) / "minigallery_test.html"
    with open(minigallery_html, "r") as fid:
        mg_html = fid.read()
    # Check thumbnail correct
    assert "_images/sphx_glr_plot_nested_thumb.png" in mg_html
    # Check href correct
    assert "sphx-glr-ex-sub-folder-sub-sub-folder-plot-nested-py" in mg_html


def _get_minigallery_thumbnails(rst_fname):
    """Check the minigallery example file number present in a rst file.

    Note this should only be used for examples in the root
    `sphinx_gallery/tests/testconfs` directory.
    """
    locator = "sphx-glr-thumbcontainer"
    regex = r".+sphx-glr-thumbcontainer.+sphx_glr_plot_(\d)_thumb.+"
    example_numbers = list()
    with open(rst_fname, "r", encoding="utf-8") as fid:
        for line in fid:
            if locator in line:
                example_numbers.append(re.match(regex, line).group(1))
    return example_numbers


@pytest.mark.add_conf(
    content="""
from sphinx_gallery.sorting import FunctionSortKey
from sphinx_gallery.utils import custom_minigallery_sort_order_sorter

sphinx_gallery_conf = {
    'examples_dirs': 'src',
    'gallery_dirs': 'ex',
    'minigallery_sort_order': FunctionSortKey(custom_minigallery_sort_order_sorter),
}"""
)
@pytest.mark.add_rst(
    file="""
Header
======

.. minigallery:: src/plot_1.py src/plot_2.py src/plot_3.py
"""
)
def test_minigallery_sort_order_callable(sphinx_app_wrapper):
    """Check `minigallery_sort_order` works when a callable."""
    sphinx_app = sphinx_app_wrapper.build_sphinx_app()

    rst_fname = Path(sphinx_app.outdir, "minigallery_test.html")
    file_numbers = _get_minigallery_thumbnails(rst_fname)
    assert file_numbers == ["3", "2", "1"]


@pytest.mark.add_conf(
    content="""
sphinx_gallery_conf = {
    'examples_dirs': 'src',
    'gallery_dirs': 'ex',
}"""
)
@pytest.mark.add_rst(
    file="""
Header
======

.. minigallery:: src/plot_1.py src/plot_1.py
"""
)
def test_minigallery_duplicate_path_input(sphinx_app_wrapper):
    """Check minigallery duplicate input paths are de-deuplicated."""
    sphinx_app = sphinx_app_wrapper.build_sphinx_app()

    rst_fname = Path(sphinx_app.outdir, "minigallery_test.html")
    file_numbers = _get_minigallery_thumbnails(rst_fname)
    assert file_numbers == ["1"]


@pytest.mark.add_conf(
    content="""
sphinx_gallery_conf = {
    'examples_dirs': 'src',
    'gallery_dirs': 'ex',
    'backreferences_dir': 'gen_modules/backreferences',
    'doc_module': ('numpy',),
}"""
)
@pytest.mark.add_rst(
    file="""
Header
======

.. minigallery:: sphinx_gallery.py_source_parser.Block src/plot_1.py
"""
)
def test_minigallery_duplicate_object_path_input(sphinx_app_wrapper):
    """Check object and path input de-deplication works in minigallery directive.

    `Block` is used in `src/plot_1.py`.
    """
    sphinx_app = sphinx_app_wrapper.build_sphinx_app()

    rst_fname = Path(sphinx_app.outdir, "minigallery_test.html")
    file_numbers = _get_minigallery_thumbnails(rst_fname)
    assert file_numbers == ["1"]


def test_write_computation_times_noop(sphinx_app_wrapper):
    app = sphinx_app_wrapper.create_sphinx_app()
    write_computation_times(app.config.sphinx_gallery_conf, None, [])


def test_write_api_usage_noop(sphinx_app_wrapper):
    write_api_entry_usage(sphinx_app_wrapper.create_sphinx_app(), list(), None)


@pytest.mark.add_conf(
    content="""
sphinx_gallery_conf = {
    'pypandoc': ['list',],
}"""
)
def test_pypandoc_config_list(sphinx_app_wrapper):
    """Tests 'pypandoc' type checking."""
    with pytest.raises(
        ConfigError,
        match=r"'pypandoc' config allowed types: \['dict', 'bool'\].",
    ):
        app = sphinx_app_wrapper.create_sphinx_app()
        fill_gallery_conf_defaults(app, app.config, check_keys=False)


@pytest.mark.add_conf(
    content="""
sphinx_gallery_conf = {
    'pypandoc': {'bad key': 1},
}"""
)
def test_pypandoc_config_keys(sphinx_app_wrapper):
    """Tests 'pypandoc' dictionary key checking."""
    with pytest.raises(
        ConfigError, match="'pypandoc' only accepts the following key values:"
    ):
        app = sphinx_app_wrapper.create_sphinx_app()
        fill_gallery_conf_defaults(app, app.config, check_keys=False)


@pytest.mark.add_conf(
    content="""
extensions += ['jupyterlite_sphinx']

sphinx_gallery_conf = {
    'backreferences_dir' : os.path.join('modules', 'gen'),
    'examples_dirs': 'src',
    'gallery_dirs': ['ex'],
}"""
)
def test_create_jupyterlite_contents(sphinx_app_wrapper):
    """Test that JupyterLite contents are created properly."""
    pytest.importorskip("jupyterlite_sphinx")
    sphinx_app = sphinx_app_wrapper.create_sphinx_app()
    gallery_conf = sphinx_app.config.sphinx_gallery_conf

    create_jupyterlite_contents(sphinx_app, exception=None)

    for i_file in ["plot_1", "plot_2", "plot_3"]:
        assert Path(
            sphinx_app.srcdir,
            "jupyterlite_contents",
            gallery_conf["gallery_dirs"][0],
            i_file + ".ipynb",
        ).exists()


@pytest.mark.add_conf(
    content="""
extensions += ['jupyterlite_sphinx']

sphinx_gallery_conf = {
    'backreferences_dir' : os.path.join('modules', 'gen'),
    'examples_dirs': 'src',
    'gallery_dirs': ['ex'],
    'jupyterlite': {'jupyterlite_contents': 'this_is_the_contents_dir'}
}"""
)
def test_create_jupyterlite_contents_non_default_contents(sphinx_app_wrapper):
    """Test that JupyterLite contents are created properly."""
    pytest.importorskip("jupyterlite_sphinx")
    sphinx_app = sphinx_app_wrapper.create_sphinx_app()
    gallery_conf = sphinx_app.config.sphinx_gallery_conf

    create_jupyterlite_contents(sphinx_app, exception=None)

    for i_file in ["plot_1", "plot_2", "plot_3"]:
        assert Path(
            sphinx_app.srcdir,
            "this_is_the_contents_dir",
            gallery_conf["gallery_dirs"][0],
            i_file + ".ipynb",
        ).exists()


@pytest.mark.add_conf(
    content="""
sphinx_gallery_conf = {
    'backreferences_dir' : os.path.join('modules', 'gen'),
    'examples_dirs': 'src',
    'gallery_dirs': ['ex'],
}"""
)
def test_create_jupyterlite_contents_without_jupyterlite_sphinx_loaded(
    sphinx_app_wrapper,
):
    """Test JupyterLite contents creation without jupyterlite_sphinx loaded."""
    pytest.importorskip("jupyterlite_sphinx")
    sphinx_app = sphinx_app_wrapper.create_sphinx_app()

    create_jupyterlite_contents(sphinx_app, exception=None)
    assert not Path(sphinx_app.srcdir, "jupyterlite_contents").exists()


@pytest.mark.add_conf(
    content="""
extensions += ['jupyterlite_sphinx']

sphinx_gallery_conf = {
    'backreferences_dir' : os.path.join('modules', 'gen'),
    'examples_dirs': 'src',
    'gallery_dirs': ['ex'],
    'jupyterlite': None,
}"""
)
def test_create_jupyterlite_contents_with_jupyterlite_disabled_via_config(
    sphinx_app_wrapper,
):
    """Test JupyterLite contents created with jupyterlite_sphinx loaded but disabled.

    JupyterLite disabled via config.
    """
    pytest.importorskip("jupyterlite_sphinx")
    sphinx_app = sphinx_app_wrapper.create_sphinx_app()

    create_jupyterlite_contents(sphinx_app, exception=None)
    assert not Path(sphinx_app.outdir, "jupyterlite_contents").exists()


@pytest.mark.add_conf(
    content="""
extensions += ['jupyterlite_sphinx']

def notebook_modification_function(notebook_content, notebook_filename):
    source = f'JupyterLite-specific change for {notebook_filename}'
    markdown_cell = {
        'cell_type': 'markdown',
        'metadata': {},
        'source': source
    }
    notebook_content['cells'] = [markdown_cell] + notebook_content['cells']


sphinx_gallery_conf = {
    'backreferences_dir' : os.path.join('modules', 'gen'),
    'examples_dirs': 'src',
    'gallery_dirs': ['ex'],
    'jupyterlite': {
        'notebook_modification_function': notebook_modification_function
    }
}"""
)
def test_create_jupyterlite_contents_with_modification(sphinx_app_wrapper):
    pytest.importorskip("jupyterlite_sphinx")
    sphinx_app = sphinx_app_wrapper.create_sphinx_app()
    gallery_conf = sphinx_app.config.sphinx_gallery_conf

    create_jupyterlite_contents(sphinx_app, exception=None)

    for i_file in ["plot_1", "plot_2", "plot_3"]:
        notebook_filename = Path(
            sphinx_app.srcdir,
            "jupyterlite_contents",
            gallery_conf["gallery_dirs"][0],
            i_file + ".ipynb",
        )
        assert notebook_filename.exists()

        with open(notebook_filename) as f:
            notebook_content = json.load(f)

        first_cell = notebook_content["cells"][0]
        assert first_cell["cell_type"] == "markdown"
        assert (
            f"JupyterLite-specific change for {notebook_filename}"
            in first_cell["source"]
        )
