# Author: Óscar Nájera
# License: 3-clause BSD
"""Testing the rst files generator."""

import os
import tempfile

import pytest
from sphinx.errors import ExtensionError

import sphinx_gallery.docs_resolv as sg


def test_embed_code_links_get_data():
    """Test that we can get data for code links."""
    sg._get_data("https://numpy.org/doc/1.18/reference")
    sg._get_data("http://scikit-learn.org/stable/")  # GZip


def test_shelve(tmpdir):
    """Test if shelve can cache and retrieve data after file is deleted."""
    test_string = "test information"
    tmp_cache = str(tmpdir)
    with tempfile.NamedTemporaryFile("w", delete=False) as f:
        f.write(test_string)
    try:
        # recovers data from temporary file and caches it in the shelve
        file_data = sg.get_data(f.name, tmp_cache)
    finally:
        os.remove(f.name)
    # tests recovered data matches
    assert file_data == test_string

    # test if cached data is available after temporary file has vanished
    assert sg.get_data(f.name, tmp_cache) == test_string


def test_parse_sphinx_docopts():
    data = """
    <script type="text/javascript">
      var DOCUMENTATION_OPTIONS = {
        URL_ROOT:    './',
        VERSION:     '2.0.2',
        COLLAPSE_INDEX: false,
        FILE_SUFFIX: '.html',
        HAS_SOURCE:  true,
        SOURCELINK_SUFFIX: '.txt'
      };
    </script>
    """
    assert sg.parse_sphinx_docopts(data) == {
        "URL_ROOT": "./",
        "VERSION": "2.0.2",
        "COLLAPSE_INDEX": False,
        "FILE_SUFFIX": ".html",
        "HAS_SOURCE": True,
        "SOURCELINK_SUFFIX": ".txt",
    }

    data_sphinx_175 = """
    <script type="text/javascript">
      var DOCUMENTATION_OPTIONS = {
        URL_ROOT: document.getElementById("documentation_options")\
                  .getAttribute('data-url_root'),
        VERSION:     '2.0.2',
        COLLAPSE_INDEX: false,
        FILE_SUFFIX: '.html',
        HAS_SOURCE:  true,
        SOURCELINK_SUFFIX: '.txt'
      };
    </script>
    """
    assert sg.parse_sphinx_docopts(data_sphinx_175) == {
        "VERSION": "2.0.2",
        "COLLAPSE_INDEX": False,
        "FILE_SUFFIX": ".html",
        "HAS_SOURCE": True,
        "SOURCELINK_SUFFIX": ".txt",
    }

    with pytest.raises(ExtensionError):
        sg.parse_sphinx_docopts("empty input")

    with pytest.raises(ExtensionError):
        sg.parse_sphinx_docopts("DOCUMENTATION_OPTIONS = ")

    with pytest.raises(ExtensionError):
        sg.parse_sphinx_docopts("DOCUMENTATION_OPTIONS = {")
