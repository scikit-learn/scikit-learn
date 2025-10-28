# Author: Chris Holdgraf
# License: 3-clause BSD
"""Testing the binder badge functionality."""

import os
import re
from copy import deepcopy
from unittest.mock import Mock

import pytest
from sphinx.application import Sphinx
from sphinx.errors import ConfigError

from sphinx_gallery.interactive_example import (
    _copy_binder_reqs,
    check_binder_conf,
    check_jupyterlite_conf,
    gen_binder_rst,
    gen_binder_url,
    gen_jupyterlite_rst,
)


def test_binder():
    """Testing binder URL generation and checks."""
    file_path = "blahblah/mydir/myfile.py"
    conf_base = {
        "binderhub_url": "http://test1.com",
        "org": "org",
        "repo": "repo",
        "branch": "branch",
        "dependencies": "../requirements.txt",
    }
    conf_base = check_binder_conf(conf_base)
    gallery_conf_base = {"gallery_dirs": "mydir", "src_dir": "blahblah"}

    url = gen_binder_url(file_path, conf_base, gallery_conf_base)
    expected = (
        "http://test1.com/v2/gh/org/repo/branch?filepath=notebooks/mydir/myfile.ipynb"
    )
    assert url == expected

    # Assert url quoted correctly
    conf0 = deepcopy(conf_base)
    special_file_path = "blahblah/mydir/files_&_stuff.py"
    conf0["branch"] = "100%_tested"
    url = gen_binder_url(special_file_path, conf0, gallery_conf_base)
    expected = (
        "http://test1.com/v2/gh/org/repo/"
        "100%25_tested?filepath=notebooks/mydir/files_%26_stuff.ipynb"
    )
    assert url == expected

    # Assert filepath prefix is added
    prefix = "my_prefix/foo"
    conf1 = deepcopy(conf_base)
    conf1["filepath_prefix"] = prefix
    url = gen_binder_url(file_path, conf1, gallery_conf_base)
    expected = (
        "http://test1.com/v2/gh/org/repo/"
        "branch?filepath={}/notebooks/"
        "mydir/myfile.ipynb"
    ).format(prefix)

    assert url == expected
    conf1.pop("filepath_prefix")

    # URL must have http
    conf2 = deepcopy(conf1)
    conf2["binderhub_url"] = "test1.com"
    with pytest.raises(ConfigError, match="did not supply a valid url"):
        url = check_binder_conf(conf2)

    # Assert missing params
    optional_keys = ["filepath_prefix", "notebooks_dir", "use_jupyter_lab"]
    for key in conf1.keys():
        if key in optional_keys:
            continue
        conf3 = deepcopy(conf1)
        conf3.pop(key)
        with pytest.raises(ConfigError, match="binder_conf is missing values"):
            url = check_binder_conf(conf3)

    # Dependencies file
    dependency_file_tests = ["requirements_not.txt", "doc-requirements.txt"]
    for ifile in dependency_file_tests:
        conf3 = deepcopy(conf1)
        conf3["dependencies"] = ifile
        with pytest.raises(
            ConfigError,
            match=r"Did not find one of `requirements.txt` " "or `environment.yml`",
        ):
            url = check_binder_conf(conf3)

    conf6 = deepcopy(conf1)
    conf6["dependencies"] = {"test": "test"}
    with pytest.raises(
        ConfigError, match="`dependencies` value should be a list of strings"
    ):
        url = check_binder_conf(conf6)

    # Missing requirements file should raise an error
    conf7 = deepcopy(conf1)
    conf7["dependencies"] = ["requirements.txt", "this_doesntexist.txt"]
    # Hack to test the case when dependencies point to a non-existing file
    # w/o needing to do a full build

    def apptmp():
        pass

    apptmp.srcdir = "/"
    with pytest.raises(ConfigError, match="Couldn't find the Binder requirements file"):
        url = _copy_binder_reqs(apptmp, conf7)

    # Check returns the correct object
    conf4 = check_binder_conf({})
    conf5 = check_binder_conf(None)
    for iconf in [conf4, conf5]:
        assert iconf == {}

    # Assert extra unknown params
    conf7 = deepcopy(conf1)
    conf7["foo"] = "blah"
    with pytest.raises(ConfigError, match="Unknown Binder config key"):
        url = check_binder_conf(conf7)

    # Assert using lab correctly changes URL
    conf_lab = deepcopy(conf_base)
    conf_lab["use_jupyter_lab"] = True
    url = gen_binder_url(file_path, conf_lab, gallery_conf_base)
    expected = (
        "http://test1.com/v2/gh/org/repo/"
        "branch?urlpath=lab/tree/notebooks/mydir/myfile.ipynb"
    )
    assert url == expected

    # Assert using static folder correctly changes URL
    conf_static = deepcopy(conf_base)
    file_path = "blahblah/mydir/myfolder/myfile.py"
    conf_static["notebooks_dir"] = "ntbk_folder"
    url = gen_binder_url(file_path, conf_static, gallery_conf_base)
    expected = (
        "http://test1.com/v2/gh/org/repo/"
        "branch?filepath=ntbk_folder/mydir/myfolder/myfile.ipynb"
    )
    assert url == expected


def test_gen_binder_rst(tmpdir):
    """Check binder rst generated correctly."""
    gallery_conf_base = {"gallery_dirs": [str(tmpdir)], "src_dir": "blahblah"}
    file_path = str(tmpdir.join("blahblah", "mydir", "myfile.py"))
    conf_base = {
        "binderhub_url": "http://test1.com",
        "org": "org",
        "repo": "repo",
        "branch": "branch",
        "dependencies": "../requirements.txt",
    }
    conf_base = check_binder_conf(conf_base)
    orig_dir = os.getcwd()
    os.chdir(str(tmpdir))
    try:
        rst = gen_binder_rst(file_path, conf_base, gallery_conf_base)
    finally:
        os.chdir(orig_dir)
    image_rst = " .. image:: images/binder_badge_logo.svg"
    target_rst = ":target: http://test1.com/v2/gh/org/repo/branch?filepath=notebooks/mydir/myfile.ipynb"  # noqa E501
    alt_rst = ":alt: Launch binder"
    assert image_rst in rst
    assert target_rst in rst
    assert alt_rst in rst
    image_fname = os.path.join(
        os.path.dirname(file_path), "images", "binder_badge_logo.svg"
    )
    assert os.path.isfile(image_fname)


@pytest.mark.parametrize("use_jupyter_lab", [True, False])
@pytest.mark.parametrize(
    "example_file",
    [
        os.path.join("example_dir", "myfile.py"),
        os.path.join("example_dir", "subdir", "myfile.py"),
    ],
)
def test_gen_jupyterlite_rst(use_jupyter_lab, example_file, tmpdir):
    """Check binder rst generated correctly."""
    gallery_conf = {
        "gallery_dirs": [str(tmpdir)],
        "src_dir": "blahblah",
        "jupyterlite": {"use_jupyter_lab": use_jupyter_lab},
    }
    file_path = str(tmpdir.join("blahblah", example_file))
    orig_dir = os.getcwd()
    os.chdir(str(tmpdir))
    try:
        rst = gen_jupyterlite_rst(file_path, gallery_conf)
    finally:
        os.chdir(orig_dir)
    image_rst = " .. image:: images/jupyterlite_badge_logo.svg"

    target_rst_template = (
        ":target: {root_url}/lite/{jupyter_part}.+index.html.+path={notebook_path}"
    )
    if "subdir" not in file_path:
        root_url = r"\.\."
        notebook_path = r"example_dir/myfile\.ipynb"
    else:
        root_url = r"\.\./\.\."
        notebook_path = r"example_dir/subdir/myfile\.ipynb"

    if use_jupyter_lab:
        jupyter_part = "lab"
    else:
        jupyter_part = "notebooks"

    target_rst = target_rst_template.format(
        root_url=root_url, jupyter_part=jupyter_part, notebook_path=notebook_path
    )
    alt_rst = ":alt: Launch JupyterLite"
    assert image_rst in rst
    assert re.search(target_rst, rst), rst
    assert alt_rst in rst
    image_fname = os.path.join(
        os.path.dirname(file_path), "images", "jupyterlite_badge_logo.svg"
    )
    assert os.path.isfile(image_fname)


def test_check_jupyterlite_conf():
    app = Mock(
        spec=Sphinx,
        config=dict(source_suffix={".rst": None}),
        extensions=[],
        srcdir="srcdir",
    )

    assert check_jupyterlite_conf(None, app) is None
    assert check_jupyterlite_conf({}, app) is None

    app.extensions = ["jupyterlite_sphinx"]
    assert check_jupyterlite_conf(None, app) is None
    assert check_jupyterlite_conf({}, app) == {
        "jupyterlite_contents": os.path.join("srcdir", "jupyterlite_contents"),
        "use_jupyter_lab": True,
        "notebook_modification_function": None,
    }

    conf = {
        "jupyterlite_contents": "this_is_the_contents_dir",
        "use_jupyter_lab": False,
    }
    expected = {
        "jupyterlite_contents": os.path.join("srcdir", "this_is_the_contents_dir"),
        "use_jupyter_lab": False,
        "notebook_modification_function": None,
    }
    assert check_jupyterlite_conf(conf, app) == expected

    def notebook_modification_function(notebook_content, notebook_filename):
        pass

    conf = {"notebook_modification_function": notebook_modification_function}

    expected = {
        "jupyterlite_contents": os.path.join("srcdir", "jupyterlite_contents"),
        "use_jupyter_lab": True,
        "notebook_modification_function": notebook_modification_function,
    }

    assert check_jupyterlite_conf(conf, app) == expected

    match = (
        "Found.+unknown keys.+another_unknown_key.+unknown_key.+"
        "Allowed keys are.+jupyterlite_contents.+"
        "notebook_modification_function.+use_jupyter_lab"
    )
    with pytest.raises(ConfigError, match=match):
        check_jupyterlite_conf(
            {"unknown_key": "value", "another_unknown_key": "another_value"}, app
        )
