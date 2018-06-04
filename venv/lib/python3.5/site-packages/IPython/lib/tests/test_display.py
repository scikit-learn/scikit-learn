"""Tests for IPython.lib.display.

"""
#-----------------------------------------------------------------------------
# Copyright (c) 2012, the IPython Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
from tempfile import NamedTemporaryFile, mkdtemp
from os.path import split, join as pjoin, dirname
import sys
try:
    import pathlib
except ImportError:
    pass

# Third-party imports
import nose.tools as nt

# Our own imports
from IPython.lib import display
from IPython.testing.decorators import skipif_not_numpy

#-----------------------------------------------------------------------------
# Classes and functions
#-----------------------------------------------------------------------------

#--------------------------
# FileLink tests
#--------------------------

def test_instantiation_FileLink():
    """FileLink: Test class can be instantiated"""
    fl = display.FileLink('example.txt')
    # TODO: remove if when only Python >= 3.6 is supported
    if sys.version_info >= (3, 6):
        fl = display.FileLink(pathlib.PurePath('example.txt'))

def test_warning_on_non_existant_path_FileLink():
    """FileLink: Calling _repr_html_ on non-existant files returns a warning
    """
    fl = display.FileLink('example.txt')
    nt.assert_true(fl._repr_html_().startswith('Path (<tt>example.txt</tt>)'))

def test_existing_path_FileLink():
    """FileLink: Calling _repr_html_ functions as expected on existing filepath
    """
    tf = NamedTemporaryFile()
    fl = display.FileLink(tf.name)
    actual = fl._repr_html_()
    expected = "<a href='%s' target='_blank'>%s</a><br>" % (tf.name,tf.name)
    nt.assert_equal(actual,expected)

def test_existing_path_FileLink_repr():
    """FileLink: Calling repr() functions as expected on existing filepath
    """
    tf = NamedTemporaryFile()
    fl = display.FileLink(tf.name)
    actual = repr(fl)
    expected = tf.name
    nt.assert_equal(actual,expected)

def test_error_on_directory_to_FileLink():
    """FileLink: Raises error when passed directory
    """
    td = mkdtemp()
    nt.assert_raises(ValueError,display.FileLink,td)

#--------------------------
# FileLinks tests
#--------------------------

def test_instantiation_FileLinks():
    """FileLinks: Test class can be instantiated
    """
    fls = display.FileLinks('example')

def test_warning_on_non_existant_path_FileLinks():
    """FileLinks: Calling _repr_html_ on non-existant files returns a warning
    """
    fls = display.FileLinks('example')
    nt.assert_true(fls._repr_html_().startswith('Path (<tt>example</tt>)'))

def test_existing_path_FileLinks():
    """FileLinks: Calling _repr_html_ functions as expected on existing dir
    """
    td = mkdtemp()
    tf1 = NamedTemporaryFile(dir=td)
    tf2 = NamedTemporaryFile(dir=td)
    fl = display.FileLinks(td)
    actual = fl._repr_html_()
    actual = actual.split('\n')
    actual.sort()
    # the links should always have forward slashes, even on windows, so replace
    # backslashes with forward slashes here
    expected = ["%s/<br>" % td,
                "&nbsp;&nbsp;<a href='%s' target='_blank'>%s</a><br>" %\
                 (tf2.name.replace("\\","/"),split(tf2.name)[1]),
                "&nbsp;&nbsp;<a href='%s' target='_blank'>%s</a><br>" %\
                 (tf1.name.replace("\\","/"),split(tf1.name)[1])]
    expected.sort()
    # We compare the sorted list of links here as that's more reliable
    nt.assert_equal(actual,expected)

def test_existing_path_FileLinks_alt_formatter():
    """FileLinks: Calling _repr_html_ functions as expected w/ an alt formatter
    """
    td = mkdtemp()
    tf1 = NamedTemporaryFile(dir=td)
    tf2 = NamedTemporaryFile(dir=td)
    def fake_formatter(dirname,fnames,included_suffixes):
        return ["hello","world"]
    fl = display.FileLinks(td,notebook_display_formatter=fake_formatter)
    actual = fl._repr_html_()
    actual = actual.split('\n')
    actual.sort()
    expected = ["hello","world"]
    expected.sort()
    # We compare the sorted list of links here as that's more reliable
    nt.assert_equal(actual,expected)

def test_existing_path_FileLinks_repr():
    """FileLinks: Calling repr() functions as expected on existing directory """
    td = mkdtemp()
    tf1 = NamedTemporaryFile(dir=td)
    tf2 = NamedTemporaryFile(dir=td)
    fl = display.FileLinks(td)
    actual = repr(fl)
    actual = actual.split('\n')
    actual.sort()
    expected = ['%s/' % td, '  %s' % split(tf1.name)[1],'  %s' % split(tf2.name)[1]]
    expected.sort()
    # We compare the sorted list of links here as that's more reliable
    nt.assert_equal(actual,expected)
    
def test_existing_path_FileLinks_repr_alt_formatter():
    """FileLinks: Calling repr() functions as expected w/ alt formatter
    """
    td = mkdtemp()
    tf1 = NamedTemporaryFile(dir=td)
    tf2 = NamedTemporaryFile(dir=td)
    def fake_formatter(dirname,fnames,included_suffixes):
        return ["hello","world"]
    fl = display.FileLinks(td,terminal_display_formatter=fake_formatter)
    actual = repr(fl)
    actual = actual.split('\n')
    actual.sort()
    expected = ["hello","world"]
    expected.sort()
    # We compare the sorted list of links here as that's more reliable
    nt.assert_equal(actual,expected)
    
def test_error_on_file_to_FileLinks():
    """FileLinks: Raises error when passed file
    """
    td = mkdtemp()
    tf1 = NamedTemporaryFile(dir=td)
    nt.assert_raises(ValueError,display.FileLinks,tf1.name)

def test_recursive_FileLinks():
    """FileLinks: Does not recurse when recursive=False
    """
    td = mkdtemp()
    tf = NamedTemporaryFile(dir=td)
    subtd = mkdtemp(dir=td)
    subtf = NamedTemporaryFile(dir=subtd)
    fl = display.FileLinks(td)
    actual = str(fl)
    actual = actual.split('\n')
    nt.assert_equal(len(actual), 4, actual)
    fl = display.FileLinks(td, recursive=False)
    actual = str(fl)
    actual = actual.split('\n')
    nt.assert_equal(len(actual), 2, actual)

@skipif_not_numpy
def test_audio_from_file():
    path = pjoin(dirname(__file__), 'test.wav')
    display.Audio(filename=path)

def test_code_from_file():
    c = display.Code(filename=__file__)
    assert c._repr_html_().startswith('<style>')
