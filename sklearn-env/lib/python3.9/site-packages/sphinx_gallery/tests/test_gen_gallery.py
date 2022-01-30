# -*- coding: utf-8 -*-
# Author: Óscar Nájera
# License: 3-clause BSD
r"""
Test Sphinx-Gallery
"""
import codecs
import os
import re
from pathlib import Path

import pytest

from sphinx.errors import ConfigError, ExtensionError, SphinxWarning
from sphinx_gallery.gen_gallery import (
    check_duplicate_filenames, check_spaces_in_filenames,
    collect_gallery_files, write_computation_times, _complete_gallery_conf)


def test_bad_config():
    """Test that bad config values are caught."""
    sphinx_gallery_conf = dict(example_dir='')
    with pytest.raises(ConfigError, match="example_dir.*did you mean 'examples_dirs'?.*"):  # noqa: E501
        _complete_gallery_conf(sphinx_gallery_conf, '', True, False)
    sphinx_gallery_conf = dict(n_subsection_order='')
    with pytest.raises(ConfigError, match=r"did you mean one of \['subsection_order', 'within_.*"):  # noqa: E501
        _complete_gallery_conf(sphinx_gallery_conf, '', True, False)


def test_default_config(sphinx_app_wrapper):
    """Test the default Sphinx-Gallery configuration is loaded

    if only the extension is added to Sphinx"""
    sphinx_app = sphinx_app_wrapper.create_sphinx_app()
    cfg = sphinx_app.config
    assert cfg.project == "Sphinx-Gallery <Tests>"
    # no duplicate values allowed The config is present already
    with pytest.raises(ExtensionError) as excinfo:
        sphinx_app.add_config_value('sphinx_gallery_conf', 'x', True)
    assert 'already present' in str(excinfo.value)


@pytest.mark.conf_file(content="""
sphinx_gallery_conf = {
    'examples_dirs': 'src',
    'gallery_dirs': 'ex',
}""")
def test_no_warning_simple_config(sphinx_app_wrapper):
    """Testing that no warning is issued with a simple config.

    The simple config only specifies input (examples_dirs) and output
    (gallery_dirs) directories.
    """
    sphinx_app = sphinx_app_wrapper.create_sphinx_app()
    cfg = sphinx_app.config
    assert cfg.project == "Sphinx-Gallery <Tests>"
    build_warn = sphinx_app._warning.getvalue()
    assert build_warn == ''


# This changed from passing the ValueError directly to
# raising "sphinx.errors.ConfigError" with "threw an exception"
@pytest.mark.parametrize('err_class, err_match', [
    pytest.param(ConfigError, 'Unknown module resetter',
                 id='Resetter unknown',
                 marks=pytest.mark.conf_file(
                     content="sphinx_gallery_conf={'reset_modules': ('f',)}")),
    pytest.param(ConfigError, 'Module resetter .* was not callab',
                 id='Resetter not callable',
                 marks=pytest.mark.conf_file(
                     content="sphinx_gallery_conf={'reset_modules': (1.,),}")),
])
def test_bad_reset(sphinx_app_wrapper, err_class, err_match):
    with pytest.raises(err_class, match=err_match):
        sphinx_app_wrapper.create_sphinx_app()


@pytest.mark.parametrize('err_class, err_match', [
    pytest.param(ConfigError, 'reset_modules_order must be a str',
                 id='Resetter unknown',
                 marks=pytest.mark.conf_file(
                     content=("sphinx_gallery_conf="
                              "{'reset_modules_order': 1,}"))),
    pytest.param(ConfigError, "reset_modules_order must be in",
                 id='reset_modules_order not valid',
                 marks=pytest.mark.conf_file(
                     content=("sphinx_gallery_conf="
                              "{'reset_modules_order': 'invalid',}"))),
])
def test_bad_reset_modules_order(sphinx_app_wrapper, err_class, err_match):
    with pytest.raises(err_class, match=err_match):
        sphinx_app_wrapper.create_sphinx_app()


@pytest.mark.parametrize('err_class, err_match', [
    pytest.param(ConfigError, 'Unknown css', id='CSS str error',
                 marks=pytest.mark.conf_file(
                     content="sphinx_gallery_conf={'css': ('foo',)}")),
    pytest.param(ConfigError, 'must be list or tuple', id="CSS type error",
                 marks=pytest.mark.conf_file(
                     content="sphinx_gallery_conf={'css': 1.}")),
])
def test_bad_css(sphinx_app_wrapper, err_class, err_match):
    with pytest.raises(err_class, match=err_match):
        sphinx_app_wrapper.create_sphinx_app()


@pytest.mark.conf_file(content="""
sphinx_gallery_conf = {
    'backreferences_dir': os.path.join('gen_modules', 'backreferences'),
    'examples_dirs': 'src',
    'gallery_dirs': 'ex',
}""")
def test_config_backreferences(sphinx_app_wrapper):
    """Test no warning is issued under the new configuration"""
    sphinx_app = sphinx_app_wrapper.create_sphinx_app()
    cfg = sphinx_app.config
    assert cfg.project == "Sphinx-Gallery <Tests>"
    assert cfg.sphinx_gallery_conf['backreferences_dir'] == os.path.join(
        'gen_modules', 'backreferences')
    build_warn = sphinx_app._warning.getvalue()
    assert build_warn == ''


def test_duplicate_files_warn(sphinx_app_wrapper):
    """Test for a warning when two files with the same filename exist."""
    sphinx_app = sphinx_app_wrapper.create_sphinx_app()

    files = ['./a/file1.py', './a/file2.py', 'a/file3.py', './b/file1.py']
    msg = ("Duplicate example file name(s) found. Having duplicate file names "
           "will break some links. List of files: {}")
    m = "['./b/file1.py']"

    # No warning because no overlapping names
    check_duplicate_filenames(files[:-1])
    build_warn = sphinx_app._warning.getvalue()
    assert build_warn == ''

    # Warning because last file is named the same
    check_duplicate_filenames(files)
    build_warn = sphinx_app._warning.getvalue()
    assert msg.format(m) in build_warn


def test_spaces_in_files_warn(sphinx_app_wrapper):
    """Test for a exception when an example filename has a space in it."""
    sphinx_app = sphinx_app_wrapper.create_sphinx_app()

    files = ['./a/file1.py', './a/file2.py', './a/file 3.py']
    msg = ("Example file name(s) with space(s) found. Having space(s) in "
           "file names will break some links. "
           "List of files: {}")
    m = "['./a/file 3.py']"

    # No warning because no filename with space
    check_spaces_in_filenames(files[:-1])
    build_warn = sphinx_app._warning.getvalue()
    assert build_warn == ''

    # Warning because last file has space
    check_spaces_in_filenames(files)
    build_warn = sphinx_app._warning.getvalue()
    assert msg.format(m) in build_warn


def _check_order(sphinx_app, key):
    index_fname = os.path.join(sphinx_app.outdir, '..', 'ex', 'index.rst')
    order = list()
    regex = '.*:%s=(.):.*' % key
    with codecs.open(index_fname, 'r', 'utf-8') as fid:
        for line in fid:
            if 'sphx-glr-thumbcontainer' in line:
                order.append(int(re.match(regex, line).group(1)))
    assert len(order) == 3
    assert order == [1, 2, 3]


@pytest.mark.conf_file(content="""
sphinx_gallery_conf = {
    'examples_dirs': 'src',
    'gallery_dirs': 'ex',
}""")
def test_example_sorting_default(sphinx_app_wrapper):
    """Test sorting of examples by default key (number of code lines)."""
    sphinx_app = sphinx_app_wrapper.create_sphinx_app()
    _check_order(sphinx_app, 'lines')


@pytest.mark.conf_file(content="""
from sphinx_gallery.sorting import FileSizeSortKey
sphinx_gallery_conf = {
    'examples_dirs': 'src',
    'gallery_dirs': 'ex',
    'within_subsection_order': FileSizeSortKey,
}""")
def test_example_sorting_filesize(sphinx_app_wrapper):
    """Test sorting of examples by filesize."""
    sphinx_app = sphinx_app_wrapper.create_sphinx_app()
    _check_order(sphinx_app, 'filesize')


@pytest.mark.conf_file(content="""
from sphinx_gallery.sorting import FileNameSortKey
sphinx_gallery_conf = {
    'examples_dirs': 'src',
    'gallery_dirs': 'ex',
    'within_subsection_order': FileNameSortKey,
}""")
def test_example_sorting_filename(sphinx_app_wrapper):
    """Test sorting of examples by filename."""
    sphinx_app = sphinx_app_wrapper.create_sphinx_app()
    _check_order(sphinx_app, 'filename')


@pytest.mark.conf_file(content="""
from sphinx_gallery.sorting import ExampleTitleSortKey
sphinx_gallery_conf = {
    'examples_dirs': 'src',
    'gallery_dirs': 'ex',
    'within_subsection_order': ExampleTitleSortKey,
}""")
def test_example_sorting_title(sphinx_app_wrapper):
    """Test sorting of examples by title."""
    sphinx_app = sphinx_app_wrapper.create_sphinx_app()
    _check_order(sphinx_app, 'title')


def test_collect_gallery_files(tmpdir, gallery_conf):
    """Test that example files are collected properly."""
    rel_filepaths = ['examples/file1.py',
                     'examples/test.rst',
                     'examples/README.txt',
                     'examples/folder1/file1.py',
                     'examples/folder1/file2.py',
                     'examples/folder2/file1.py',
                     'tutorials/folder1/subfolder/file1.py',
                     'tutorials/folder2/subfolder/subsubfolder/file1.py']

    abs_paths = [tmpdir.join(rp) for rp in rel_filepaths]
    for ap in abs_paths:
        ap.ensure()

    examples_path = tmpdir.join('examples')
    dirs = [examples_path.strpath]
    collected_files = set(collect_gallery_files(dirs, gallery_conf))
    expected_files = set(
        [ap.strpath for ap in abs_paths
         if re.search(r'examples.*\.py$', ap.strpath)])

    assert collected_files == expected_files

    tutorials_path = tmpdir.join('tutorials')
    dirs = [examples_path.strpath, tutorials_path.strpath]
    collected_files = set(collect_gallery_files(dirs, gallery_conf))
    expected_files = set(
        [ap.strpath for ap in abs_paths if re.search(r'.*\.py$', ap.strpath)])

    assert collected_files == expected_files


def test_collect_gallery_files_ignore_pattern(tmpdir, gallery_conf):
    """Test that ignore pattern example files are not collected."""
    rel_filepaths = ['examples/file1.py',
                     'examples/folder1/fileone.py',
                     'examples/folder1/file2.py',
                     'examples/folder2/fileone.py']

    abs_paths = [tmpdir.join(rp) for rp in rel_filepaths]
    for ap in abs_paths:
        ap.ensure()

    gallery_conf['ignore_pattern'] = r'one'
    examples_path = tmpdir.join('examples')
    dirs = [examples_path.strpath]
    collected_files = set(collect_gallery_files(dirs, gallery_conf))
    expected_files = set(
        [ap.strpath for ap in abs_paths
         if re.search(r'one', os.path.basename(ap.strpath)) is None])

    assert collected_files == expected_files


@pytest.mark.conf_file(content="""
sphinx_gallery_conf = {
    'backreferences_dir' : os.path.join('modules', 'gen'),
    'examples_dirs': 'src',
    'gallery_dirs': ['ex'],
    'binder': {'binderhub_url': 'http://test1.com', 'org': 'org',
               'repo': 'repo', 'branch': 'branch',
               'notebooks_dir': 'ntbk_folder',
               'dependencies': 'requirements.txt'}
}""")
def test_binder_copy_files(sphinx_app_wrapper, tmpdir):
    """Test that notebooks are copied properly."""
    from sphinx_gallery.binder import copy_binder_files
    sphinx_app = sphinx_app_wrapper.create_sphinx_app()
    gallery_conf = sphinx_app.config.sphinx_gallery_conf
    # Create requirements file
    with open(os.path.join(sphinx_app.srcdir, 'requirements.txt'), 'w'):
        pass
    copy_binder_files(sphinx_app, None)

    for i_file in ['plot_1', 'plot_2', 'plot_3']:
        assert os.path.exists(os.path.join(
            sphinx_app.outdir, 'ntbk_folder', gallery_conf['gallery_dirs'][0],
            i_file + '.ipynb'))


@pytest.mark.conf_file(content="""
sphinx_gallery_conf = {
    'examples_dirs': 'src',
    'gallery_dirs': 'ex',
}""")
def test_failing_examples_raise_exception(sphinx_app_wrapper):
    example_dir = os.path.join(sphinx_app_wrapper.srcdir,
                               'src')
    with codecs.open(os.path.join(example_dir, 'plot_3.py'), 'a',
                     encoding='utf-8') as fid:
        fid.write('raise SyntaxError')
    with pytest.raises(ExtensionError) as excinfo:
        sphinx_app_wrapper.build_sphinx_app()
    assert "Unexpected failing examples" in str(excinfo.value)


@pytest.mark.conf_file(content="""
sphinx_gallery_conf = {
    'examples_dirs': 'src',
    'gallery_dirs': 'ex',
    'filename_pattern': 'plot_1.py',
}""")
def test_expected_failing_examples_were_executed(sphinx_app_wrapper):
    """Testing that no exception is issued when broken example is not built

    See #335 for more details.
    """
    sphinx_app_wrapper.build_sphinx_app()


@pytest.mark.conf_file(content="""
sphinx_gallery_conf = {
    'examples_dirs': 'src',
    'gallery_dirs': 'ex',
    'only_warn_on_example_error': True,
}""")
def test_only_warn_on_example_error(sphinx_app_wrapper):
    """
    Test behaviour of only_warn_on_example_error flag.
    """
    example_dir = Path(sphinx_app_wrapper.srcdir) / 'src'
    with codecs.open(example_dir / 'plot_3.py', 'a', encoding='utf-8') as fid:
        fid.write('raise ValueError')
    sphinx_app = sphinx_app_wrapper.build_sphinx_app()

    build_warn = sphinx_app._warning.getvalue()
    assert 'plot_3.py failed to execute correctly' in build_warn
    assert 'WARNING: Here is a summary of the problems' in build_warn


@pytest.mark.conf_file(content="""
sphinx_gallery_conf = {
    'examples_dirs': 'src',
    'gallery_dirs': 'ex',
    'only_warn_on_example_error': True,
}""")
def test_only_warn_on_example_error_sphinx_warning(sphinx_app_wrapper):
    """
    Test behaviour of only_warn_on_example_error flag.
    """
    sphinx_app_wrapper.kwargs['warningiserror'] = True
    example_dir = Path(sphinx_app_wrapper.srcdir) / 'src'
    with codecs.open(example_dir / 'plot_3.py', 'a', encoding='utf-8') as fid:
        fid.write('raise ValueError')
    with pytest.raises(SphinxWarning) as excinfo:
        sphinx_app_wrapper.build_sphinx_app()
    assert "plot_3.py failed to execute" in str(excinfo.value)


@pytest.mark.conf_file(content="""
sphinx_gallery_conf = {
    'examples_dirs': 'src',
    'gallery_dirs': 'ex',
    'expected_failing_examples' :['src/plot_2.py'],
}""")
def test_examples_not_expected_to_pass(sphinx_app_wrapper):
    with pytest.raises(ExtensionError) as excinfo:
        sphinx_app_wrapper.build_sphinx_app()
    assert "expected to fail, but not failing" in str(excinfo.value)


@pytest.mark.conf_file(content="""
sphinx_gallery_conf = {
    'show_memory': lambda func: (0., func()),
    'gallery_dirs': 'ex',
}""")
def test_show_memory_callable(sphinx_app_wrapper):
    sphinx_app = sphinx_app_wrapper.build_sphinx_app()
    status = sphinx_app._status.getvalue()
    assert "0.0 MB" in status


@pytest.mark.conf_file(content="""
sphinx_gallery_conf = {
    'first_notebook_cell': 2,
}""")
def test_first_notebook_cell_config(sphinx_app_wrapper):
    from sphinx_gallery.gen_gallery import parse_config
    # First cell must be str
    with pytest.raises(ConfigError):
        parse_config(sphinx_app_wrapper.create_sphinx_app(), False)


@pytest.mark.conf_file(content="""
sphinx_gallery_conf = {
    'last_notebook_cell': 2,
}""")
def test_last_notebook_cell_config(sphinx_app_wrapper):
    from sphinx_gallery.gen_gallery import parse_config
    # First cell must be str
    with pytest.raises(ConfigError):
        parse_config(sphinx_app_wrapper.create_sphinx_app(), False)


@pytest.mark.conf_file(content="""
sphinx_gallery_conf = {
    'backreferences_dir': False,
}""")
def test_backreferences_dir_config(sphinx_app_wrapper):
    """Tests 'backreferences_dir' type checking."""
    from sphinx_gallery.gen_gallery import parse_config
    with pytest.raises(ConfigError,
                       match="The 'backreferences_dir' parameter must be of"):
        parse_config(sphinx_app_wrapper.create_sphinx_app(), False)


@pytest.mark.conf_file(content="""
import pathlib

sphinx_gallery_conf = {
    'backreferences_dir': pathlib.Path('.'),
}""")
def test_backreferences_dir_pathlib_config(sphinx_app_wrapper):
    """Tests pathlib.Path does not raise exception."""
    from sphinx_gallery.gen_gallery import parse_config
    parse_config(sphinx_app_wrapper.create_sphinx_app(), False)


def test_write_computation_times_noop():
    write_computation_times(None, None, [[[0]]])


@pytest.mark.conf_file(content="""
sphinx_gallery_conf = {
    'pypandoc': ['list',],
}""")
def test_pypandoc_config_list(sphinx_app_wrapper):
    """Tests 'pypandoc' type checking."""
    from sphinx_gallery.gen_gallery import parse_config
    with pytest.raises(ConfigError,
                       match="'pypandoc' parameter must be of type bool or "
                             "dict"):
        parse_config(sphinx_app_wrapper.create_sphinx_app(), False)


@pytest.mark.conf_file(content="""
sphinx_gallery_conf = {
    'pypandoc': {'bad key': 1},
}""")
def test_pypandoc_config_keys(sphinx_app_wrapper):
    """Tests 'pypandoc' dictonary key checking."""
    from sphinx_gallery.gen_gallery import parse_config
    with pytest.raises(ConfigError,
                       match="'pypandoc' only accepts the following key "
                             "values:"):
        parse_config(sphinx_app_wrapper.create_sphinx_app(), False)
